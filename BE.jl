#—————————————————————————————————————————————CALIBRATION——————————————————————————————————————————
using GEMS
import GEMS.transmission_probability
using Parameters
using Distributions 
using DataFrames
using Plots

@with_kw mutable struct SettingRate <:GEMS.TransmissionFunction   
    Household_rate::Float64
    OfficeSchool_rate::Float64
    General_rate::Float64
end

function GEMS.transmission_probability(transFunc::SettingRate, infecter::Individual, infected::Individual, setting::Setting, tick::Int16)
    if isa(setting, Household)
        return transFunc.Household_rate
    elseif isa(setting, Workplace) || isa(setting, WorkplaceSite) || isa(setting, Office) || isa(setting, Department) || isa(setting, SchoolComplex) || isa(setting, SchoolYear) || isa(setting, School) || isa(setting, SchoolClass)
        return transFunc.OfficeSchool_rate
    else
        return transFunc.General_rate
    end
end

dp_h = Hospitalized(exposure_to_infectiousness_onset = Poisson(3), 					           infectiousness_onset_to_symptom_onset = Poisson(1),            		               	symptom_onset_to_severeness_onset = Poisson(2),            			           severeness_onset_to_hospital_admission = Poisson(1),            	              		hospital_admission_to_hospital_discharge = Poisson(10),            	              		hospital_discharge_to_severeness_offset = Poisson(2),            	            		severeness_offset_to_recovery = Poisson(1)); 

dp_s = Symptomatic(exposure_to_infectiousness_onset = Poisson(3),
			infectiousness_onset_to_symptom_onset = Poisson(1),
			symptom_onset_to_recovery = Poisson(7));    

aba = AgeBasedProgressionAssignment(age_groups = ["-59", "60-"], 
						progression_categories = ["Symptomatic", "Hospitalized"],        						stratification_matrix=[[0.1, 0.9],
									  [0.3, 0.7]]);  
                                      
tf = SettingRate(Household_rate=0.4, OfficeSchool_rate = 0.125, General_rate = 0.32)

p = Pathogen(id = 1,             
		name = "Covid19",             
		progressions = [dp_s, dp_h],             
		progression_assignment = aba,             
		transmission_function = tf);  


#function to calculate share of each category in total infections
function setting_shares(inf_df)
           h = sum(inf_df.setting_type .== 'h')
           g = sum(inf_df.setting_type .== '?' .|| inf_df.setting_type .== 'm')
	o = sum(inf_df.setting_type .!= 'h' .&& inf_df.setting_type .!= '?' .&& inf_df.setting_type .!= 'm')

           total = h + o + g

           return Dict(
               :Household => h/total,
               :OfficeSchool => o/total,
               :General => g/total
           )
end

#function to calculate share of people older than 60 in total hospitalization
function infection_share_60plus(inf_df, ϵ = 1e-6)
    num = sum(inf_df.hospital_admission .!= -1 .&& inf_df.age_b .> 60)
    den = sum(inf_df.hospital_admission .!= -1)
    return num / max(den, ϵ)
end

#—————————————————————————BE—————————————————————

tf_be = SettingRate(Household_rate=0.25, OfficeSchool_rate = 0.11, General_rate = 0.24);
aba_be = AgeBasedProgressionAssignment(age_groups = ["-59", "60-"], 
						progression_categories = ["Symptomatic", "Hospitalized"],        						
                        stratification_matrix=[[0.8, 0.2],
									            [0.01, 0.99]]);
aba_be.stratification_matrix = [[0.8, 0.2] [0.01, 0.99]]
p_be = Pathogen(id = 1,             
		name = "Covid19",             
		progressions = [dp_s, dp_h],             
		progression_assignment = aba_be,             
		transmission_function = tf_be);  

custom_be = Simulation(population = "BE", pathogen = p_be, "label = without intervention");
run!(custom_be);

inf_df_be = infectionsDF(PostProcessor(custom_be));
obs_be = setting_shares(inf_df_be)
obs60_be = infection_share_60plus(inf_df_be)

println(obs_be)
println(obs60_be)

#——————————————————————————————————HALVING TOTAL INFECTIONS———————————————————————————————————————

# function to calculate ratio between total infections with and without intervention
function calc_inf_ratio(inf_df, sim_wi)
	#inf_df = infectionsDF(PostProcessor(sim))
	old_count = size(inf_df)[1]
	inf_df_wi = infectionsDF(PostProcessor(sim_wi))
	current_count = size(inf_df_wi)[1]
	return(current_count/old_count)
end

#————————————————————BE———————————————————————----------------------------------------------

be_inf_wi = Simulation(population = "BE", pathogen = p_be, label = "Berlin with intervention(inf)")
self_isolation = IStrategy("Self Isolation", be_inf_wi)
add_measure!(self_isolation, SelfIsolation(14))
# identify list of household members
# find_household_members = IStrategy("Find Household Members", be_inf_wi)
# add_measure!(find_household_members, FindSettingMembers(Household, self_isolation))
# trigger "household tracing"
# trigger = SymptomTrigger(find_household_members)
trigger = SymptomTrigger(self_isolation)
add_symptom_trigger!(be_inf_wi, trigger)

school_closing = SStrategy("Close Schools", be_inf_wi)
add_measure!(school_closing, CloseSetting())

trigger = STickTrigger(School, school_closing, switch_tick = Int16(10))
add_tick_trigger!(be_inf_wi, trigger)

office_closing = SStrategy("close-offices", be_inf_wi)
add_measure!(office_closing, CloseSetting())

trigger = STickTrigger(Office, office_closing, switch_tick = Int16(5))
add_tick_trigger!(be_inf_wi, trigger)


run!(be_inf_wi);

inf_ratio_be = calc_inf_ratio(inf_df_be, be_inf_wi)
println(inf_ratio_be)

rd_be_woi = ResultData(custom_be)
rd_be_inf_wi = ResultData(be_inf_wi)

p_inf_be = gemsplot([rd_be_woi, rd_be_inf_wi], type = (:CumulativeCases, :TickCases))

#——————————————————————————————————HALVING TOTAL HOSPITALIZATION——————————————————————————————————————
# function to calculate ratio between total hospitalizations with and without intervention
function calc_hsp_ratio(sim, sim_wi)
	inf_df = infectionsDF(PostProcessor(sim))
	old_count = sum(inf_df.hospital_admission .!= -1)
	inf_df_wi = infectionsDF(PostProcessor(sim_wi))
	current_count = sum(inf_df_wi.hospital_admission .!= -1)
	return(current_count/old_count)
end

#————————————————————————BE——————————————————————————————
be_h_wi = Simulation(population = "BE", pathogen = p_be, label = "Berlin with intervention(h)")
self_isolation = IStrategy("Self Isolation", be_h_wi)
add_measure!(self_isolation, SelfIsolation(14))
# identify list of household members
#find_household_members = IStrategy("Find Household Members", be_h_wi)
#add_measure!(find_household_members, FindSettingMembers(Household, self_isolation))
# trigger "household tracing"
trigger = SymptomTrigger(self_isolation)
#trigger = SymptomTrigger(find_household_members)
add_symptom_trigger!(be_h_wi, trigger)

reduce_contacts = SStrategy("Reduce Contacts", be_h_wi)
add_measure!(reduce_contacts, ChangeContactMethod(ContactparameterSampling(0.7)))

trigger = STickTrigger(SchoolClass, reduce_contacts, switch_tick = Int16(1))
add_tick_trigger!(be_h_wi, trigger)

school_closing = SStrategy("Close Schools", be_h_wi)
add_measure!(school_closing, CloseSetting())

trigger = STickTrigger(School, school_closing, switch_tick = Int16(1))
add_tick_trigger!(be_h_wi, trigger)

office_closing = SStrategy("close-offices", be_h_wi)
add_measure!(office_closing, CloseSetting())

trigger = STickTrigger(Office, office_closing, switch_tick = Int16(1))
add_tick_trigger!(be_h_wi, trigger)

run!(be_h_wi);

hsp_ratio_be = calc_hsp_ratio(custom_be, be_h_wi);
println(hsp_ratio_be)

#rd_be_woi = ResultData(custom_be)
rd_be_h_wi = ResultData(be_h_wi)

p_be_h = gemsplot([rd_be_woi, rd_be_h_wi], type = (:HospitalOccupancy))