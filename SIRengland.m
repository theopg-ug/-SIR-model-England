%% SIR deterministic model simulations for england
T = readtable("nation_newAdmissions.csv"); % specifies the file 

England = strcmp(strtrim(T.area_name), 'England'); % removes white space to find all the rows containing england
T(England, :); % specifies all columns
England_data = T(England, :); % stores england data so it can be plotted
England_dates = datetime(England_data.date, 'InputFormat', 'yyyy-MM-dd'); % changes how the x-axis is plotted

pop_size = 56e6; % number of people in England
I = 1; % initial number of infected and infectious in the population
R = 0; % initial number of recovered and immune in the population
S = pop_size - I - R; % initial number of susceptibles in the population (everyone else is susceptible)
initial_conditions = [S,I,R]/pop_size; % put together the initial conditions as a proportion of the population

params.beta=0.275; % transmission rate
params.gamma=0.1; % per day recovery rate as stated in the literature
start_date = datetime(2020,2,15); %Specifies when the SIR model should start

time_interval = 0:1:101; % output times, outputs once a day for 101 days
[t,y] = ode45(@(t,y) SIR_model(t,y,params),time_interval,initial_conditions); % run the simulation, using the function below. ode45 takes values at a specific time from dydt, they get plotted and then moves forward in time by 1 and passes to SIR to find the value for that time
model_dates = start_date + days(t) ; % converts the days into real dates from the data but cannot happen until t has been created by ode45
lockdown = datetime(2020,4,1); % sets lockdown date at the infection peak, moved to stop initial curve being exponential 
t_lockdown = days(lockdown - start_date) + 1; % calculates time of lockdown but need to add one as this only calculates dates in between
figure('Name','England') % Creates a figure name, making it easier to distinguish when lots of figures are open
plot(model_dates(1:t_lockdown),y(1:t_lockdown,2)*pop_size,"Color","b") % create an empty figure and plot I over the model dates
title('Hospital Admissions during the first Covid-19 Wave in England') % Adds the title to the figure
hold on % keeps the SIR model plot and adds whatever is plotted next
plot(England_dates, England_data.value) % plots the england hospital admissions data on top of the SIR plot
ylim([0,3500]) % sets a Y-axis limit to make the data easier to digest
xlim([datetime(2020,1,1),datetime(2020,7,1)]) % sets an x-axis limit to make the data easier to digest
xlabel('Time (months)'), ylabel('Hospital Admissions (per day)') % label axes

S_start = y(1,1); % stores the value for S at the start of the simulation so total infections be found
params.beta = 0.064; % fits the secondary model (beta 1) line to the data
new_init = y(model_dates==lockdown,:); % Sets a new initial conditions where the SIR values from the given lockdown date are taken and provided as the starting SIR values
[t,y] = ode45(@(t,y) SIR_model(t,y,params),time_interval,new_init); % run the simulation, using the function below. ode45 takes values at a specific time from dydt, they get plotted and then moves forward in time by 1 and passes to SIR to find the value for that time
new_model_dates = lockdown + days(t); % starts the lockdown line from the date stated and follows on from there for the duration stated
plot(new_model_dates,y(:,2)*pop_size, "Color","b") % plots the declining line based on the infectious population relative to the nations population in a specific colour
legend ('SIR model','England Hospitalisations') % labels the data points
hold off % nothing else to plot 

%%%%%%%% This function gives the SIR model equations
function dydt = SIR_model(t,y,params) % the function that ode45 is using data from
S = y(1); % defines Susceptible population for use in future equations
I = y(2); % defines Infectious population for use in future equations
R = y(3); % defines Recovered population for use in future equations

beta = params.beta; % creates variables to use relative to the previous variables
gamma = params.gamma; % creates variables to use relative to the previous variables

dSdt = -beta*S*I; % rate of change of S
dIdt = beta*S*I - gamma*I; % Rate of change from S to I and I to R
dRdt = gamma*I; % Rate of change from I to R

dydt = [dSdt;dIdt;dRdt]; % recall the values into y and t so the ode45 can step forward in time 

end