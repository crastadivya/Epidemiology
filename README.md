# Epidemiology
Simulated infections data for Berlin and Saxony-Anhalt states using GEMS in Julia such that certain assumptions are satisfied. Then interventions are suggested to makethe infections and hospitalizations into half.

The following assumptions are considered:

## Assumptions related to infections:
  * 1/3 of infections come from Household
  * 1/3 of infections come from school and work place
  * 1/3 of infections come from other places

## Assumptions related to hospitalizations:
  * 50½ of hospitalizations comes from people aged > 60

## Goal of interventions:
  * Reduce the total infections to half
  * Reduce the total hospitalizations to half

## Description of files:

* BE.jl and ST.jl contains code written in Julia to simulate the data with above assumptions and implement intervention to achieve goals for Berlin and Saxony-Anhalt respectively. 
* BE_inf and BE_hsp shows the results of infections and hospitalizations with and without intervention for Berlin.
* ST_inf and ST_hsp shows the results of infections and hospitalizations with and without intervention for Saxony Anhalt.
* Report.pdf is the report with problem description, methods, result, and analysis.
* Presentation.pdf is the PPT presentation of Report.
