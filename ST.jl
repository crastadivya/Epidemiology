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

#—————————————————————————ST—————————————————————

tf_st = SettingRate(Household_rate=0.4, OfficeSchool_rate = 0.125, General_rate = 0.32)
aba_st = AgeBasedProgressionAssignment(age_groups = ["-59", "60-"], 
						progression_categories = ["Symptomatic", "Hospitalized"],        						
                        stratification_matrix=[[0.75, 0.25],
									            [0.03, 0.97]]);

aba_st.stratification_matrix = [[0.75, 0.25] [0.03, 0.97]]
p_st = Pathogen(id = 1,             
		name = "Covid19",             
		progressions = [dp_s, dp_h],             
		progression_assignment = aba_st,             
		transmission_function = tf_st);  

custom_st = Simulation(population = "ST",  pathogen = p_st, label = "ST without intervention");
run!(custom_st);

inf_df_st = infectionsDF(PostProcessor(custom_st));
obs_st = setting_shares(inf_df_st)
obs60_st = infection_share_60plus(inf_df_st)

println(obs_st)
println(obs60_st)

#——————————————————————————————————HALVING TOTAL INFECTION——————————————————————————————————————
# function to calculate ratio between total infections with and without intervention
function calc_inf_ratio(sim, sim_wi)
	inf_df = infectionsDF(PostProcessor(sim))
	old_count = size(inf_df)[1]
	inf_df_wi = infectionsDF(PostProcessor(sim_wi))
	current_count = size(inf_df_wi)[1]
	return(current_count/old_count)
end

#————————————————————ST———————————————————————-------------------------------------------

st_inf_wi = Simulation(population = "ST", pathogen = p_st, label = "ST with intervention(inf)")
self_isolation = IStrategy("Self Isolation", st_inf_wi)
add_measure!(self_isolation, SelfIsolation(7))
# identify list of household members
find_household_members = IStrategy("Find Household Members", st_inf_wi)
add_measure!(find_household_members, FindSettingMembers(Household, self_isolation))
# trigger "household tracing"
trigger = SymptomTrigger(find_household_members)
add_symptom_trigger!(st_inf_wi, trigger)

office_closing = SStrategy("close-offices", st_inf_wi)
add_measure!(office_closing, CloseSetting())

trigger = STickTrigger(Office, office_closing, switch_tick = Int16(10))
add_tick_trigger!(st_inf_wi, trigger)

run!(st_inf_wi);

inf_ratio_st = calc_inf_ratio(custom_st, st_inf_wi)
println(inf_ratio_st)

rd_st_woi = ResultData(custom_st)
rd_st_inf_wi = ResultData(st_inf_wi)

p_inf_st = gemsplot([rd_st_woi, rd_st_inf_wi], type = (:CumulativeCases, :TickCases), legend=:topleft)


#——————————————————————————————————HALVING TOTAL HOSPITALIZATION——————————————————————————————————————
# function to calculate ratio between total hospitalizations with and without intervention
function calc_hsp_ratio(sim, sim_wi)
	inf_df = infectionsDF(PostProcessor(sim))
	old_count = sum(inf_df.hospital_admission .!= -1)
	inf_df_wi = infectionsDF(PostProcessor(sim_wi))
	current_count = sum(inf_df_wi.hospital_admission .!= -1)
	return(current_count/old_count)
end


#————————————————————————ST——————————————————————————————
st_h_wi = Simulation(population = "ST", pathogen = p_st, label = "ST with intervention(h)")
self_isolation = IStrategy("Self Isolation", st_h_wi)
add_measure!(self_isolation, SelfIsolation(7))
# identify list of household members
find_household_members = IStrategy("Find Household Members", st_h_wi)
add_measure!(find_household_members, FindSettingMembers(Household, self_isolation))
# trigger "household tracing"
trigger = SymptomTrigger(find_household_members)
add_symptom_trigger!(st_h_wi, trigger)

reduce_contacts = SStrategy("Reduce Contacts", st_h_wi)
add_measure!(reduce_contacts, ChangeContactMethod(ContactparameterSampling(0.25)))

trigger = STickTrigger(SchoolClass, reduce_contacts, switch_tick = Int16(100))
add_tick_trigger!(st_h_wi, trigger)

trigger = STickTrigger(Office, reduce_contacts, switch_tick = Int16(100))
add_tick_trigger!(st_h_wi, trigger)

school_closing = SStrategy("Close Schools", st_h_wi)
add_measure!(school_closing, CloseSetting())

trigger = STickTrigger(School, school_closing, switch_tick = Int16(100))
add_tick_trigger!(st_h_wi, trigger)

office_closing = SStrategy("close-offices", st_h_wi)
add_measure!(office_closing, CloseSetting())

trigger = STickTrigger(Office, office_closing, switch_tick = Int16(50))
add_tick_trigger!(st_h_wi, trigger)


run!(st_h_wi);

hsp_ratio_st = calc_hsp_ratio(custom_st, st_h_wi);
println(hsp_ratio_st)

#rd_st_woi = ResultData(custom_st)
rd_st_h_wi = ResultData(st_h_wi)

p_h_st = gemsplot([rd_st_woi, rd_st_h_wi], type = (:HospitalOccupancy), legend=:outerbottom)