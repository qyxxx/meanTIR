import numpy as np
import pandas as pd
from lifelines import CoxTimeVaryingFitter

def naive_est(data, min_time = 0, max_time = 1440*7-5,boot = 'NULL', id_col = 'patient_id', time = 'time', value_in_range = 'value_in_range'):
  #calculate TIR within max time selected
  data = data.loc[data[time] <= max_time]
  data = data.loc[data[time] >= min_time]
  #naive method groupby patient id calcualte TIR_i then average
  TIR = (data.groupby(id_col)[value_in_range].mean().reset_index())[value_in_range].mean()
  if not(str(boot).isnumeric()):
    return TIR
  else:
    boot_TIR = []
    for i in range(boot):
      boot_x = data.merge(pd.DataFrame(np.random.choice(list(data[id_col].unique()), size=len(list(data[id_col].unique())), replace=True), columns = [id_col]))
      boot_x['repeat_id_num'] = boot_x.groupby([id_col, time]).cumcount().add(1).apply(str)
      boot_x[id_col] = boot_x[id_col] + '_'+ boot_x['repeat_id_num']
      boot_TIR.append(naive_est(boot_x, min_time = min_time, max_time = max_time, boot = 'NULL', id_col=id_col, time=time, value_in_range=value_in_range))
    return [TIR, np.std(boot_TIR, ddof=1)]
  
def proposed_est_noninfo(data, min_time = 0, max_time = 1440*7-5, boot = 'NULL', id_col = 'patient_id', time = 'time', value_in_range = 'value_in_range'):
  data = data.loc[data[time] <= max_time]
  data = data.loc[data[time] >= min_time]
  TIR = (data.groupby(time)[value_in_range].mean().reset_index())[value_in_range].mean()
  if not(str(boot).isnumeric()):
    return TIR
  else:
    boot_TIR = []
    for i in range(boot):
      boot_x = data.merge(pd.DataFrame(np.random.choice(list(data[id_col].unique()), size=len(list(data[id_col].unique())), replace=True), columns = [id_col]))
      boot_x['repeat_id_num'] = boot_x.groupby([id_col, time]).cumcount().add(1).apply(str)
      boot_x[id_col] = boot_x[id_col] + '_'+ boot_x['repeat_id_num']
      boot_TIR.append(proposed_est_noninfo(boot_x, min_time = min_time, max_time = max_time, boot = 'NULL', id_col=id_col, time=time, value_in_range=value_in_range))
    return [TIR, np.std(boot_TIR, ddof=1)]

def proposed_est_cox(data, min_time = 0, max_time = 1440*7-5, id_col="patient_id", event_col="event", start_col="time", stop_col="time2", formula='var1', boot = 'NULL', value_in_range = 'value_in_range'):
  dat = data
  time = start_col
  patient_id = id_col
  #fit cox model with time dep covariate
  ctv = CoxTimeVaryingFitter(penalizer=0.1)
  ctv.fit(dat, id_col = id_col, event_col= event_col, start_col=start_col, stop_col=stop_col, formula=formula)
  cumulative_hazard = ctv.baseline_cumulative_hazard_
  cumulative_hazard[time] = cumulative_hazard.index
  cumulative_hazard['hazard diff'] = cumulative_hazard['baseline hazard'].diff().fillna(cumulative_hazard['baseline hazard'])
  dat = dat.reset_index(drop=True)
  dat['predict partial hazard'] = ctv.predict_partial_hazard(dat).reset_index(drop=True)
  dat = dat.merge(cumulative_hazard, on = time, how = 'left').fillna(value = 0)
  dat['lambda_exp_diff'] = dat['hazard diff'] * dat['predict partial hazard']
  dat['cum_lambda_exp_diff'] = dat.groupby([patient_id])['lambda_exp_diff'].cumsum()
  dat['weight'] = 1/(np.exp(- dat['cum_lambda_exp_diff']))
  #calculate TIR within max time selected
  dat = dat.loc[dat[time] <= max_time]
  dat = dat.loc[dat[time] >= min_time]
  TIR = (dat.groupby(time))[[value_in_range, 'weight']].apply(lambda x: np.average(x[value_in_range], weights = x['weight'])).mean()
  #boots or not
  if not(str(boot).isnumeric()):
    return TIR
  else:
    boot_TIR = []
    for i in range(boot):
      boot_x = data.merge(pd.DataFrame(np.random.choice(list(data[patient_id].unique()), size=len(list(data[patient_id].unique())), replace=True), columns = [patient_id]))
      boot_x['repeat_id_num'] = boot_x.groupby([patient_id, time]).cumcount().add(1).apply(str)
      boot_x[patient_id] = boot_x[patient_id] + '_'+ boot_x['repeat_id_num']
      boot_TIR.append(proposed_est_cox(boot_x, min_time = min_time, max_time = max_time, id_col = id_col, event_col= event_col, start_col=start_col, stop_col=stop_col, formula=formula, boot = 'NULL', value_in_range=value_in_range))
    return [TIR, np.std(boot_TIR, ddof=1)]