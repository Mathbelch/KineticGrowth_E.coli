import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys as sys

## Read data from a txt file in csv formmat and return an np.array with time and a Media Object, where keys are medium number and values are np.arrays with doptical density for each medium:
def getting_object(file_path):
    
   time = []
   Media = {}

   with open(file_path, encoding="utf-8") as f:
       conteudo = f.readlines() 
       f.close
       dados = []
   for line in range(len(conteudo)):
       dados.append(conteudo[line].split(','))
   for row in range(1, len(dados[0])):
       Media[f"{row-1}"] = []
   for line in range(len(dados)):
       time.append(float(dados[line][0]))
       for row in range(1, len(dados[line])):
           Media[str(row-1)].append(float(dados[line][row]))
   for medium in range(len(Media)):
      Media[str(medium)] = np.array(Media[str(medium)])
   time_array = np.array(time)
   return (time_array, Media)
  
## Checks DO values for maximum especific rate (muMax) and compare other values with it whithin a range stablished. Consecutives values inside this range are considered to be part of the exponential phase of growth.
def finding_exponential_growth(time_interval, Media, time, faixa):
   
   exponential_growth = {}
   
   for medium in range(len(Media)):
      DO = np.array(Media[str(medium)])
      DO = np.log10(DO)
      DO = np.diff(DO)
      instantMu = ((DO*2.303)/time_interval)
      muMax = np.max(instantMu)
      time_muMax_Index = np.where(instantMu == muMax)[0][0]
      exponention_growth_indexes = []
      
      for i in range(time_muMax_Index - 1, -1, -1):
         if(instantMu[i] > (muMax * (1 - faixa[medium]))):
            exponention_growth_indexes.insert(0,i)
         else:
            break
         
      exponention_growth_indexes.append(time_muMax_Index)
      exponention_growth_indexes.append(time_muMax_Index+1)
      
      for i in range(time_muMax_Index + 1, len(instantMu)-2):
         if(instantMu[i] > (muMax * (1 - faixa[medium]))):
            exponention_growth_indexes.append(i+1)
         else: 
            break
         
      if len(exponention_growth_indexes) < 3:
         print(f'São necessários ao menos três pontos para o ajuste das curvas do meio {medium}, tenta aumentar a faixa de variação da velocidade específica, que atualmente está em {faixa[medium]}.')
         sys.exit()
         
      exponential_growth[str(medium)] = np.array(exponention_growth_indexes)
      
   return exponential_growth

## Obtaining Kinetic parameters base on Matematical approach.
def obtaining_kinetic_parameters_Matematical(Media, time, exponential_growth):
   
   results = {}
   
   for medium in range(len(Media)):
      DO1 = np.log10(Media[str(medium)][exponential_growth[str(medium)][0]])
      DO2 = np.log10(Media[str(medium)][exponential_growth[str(medium)][-1]])
      time_exp = [time[exponential_growth[str(medium)][0]], time[exponential_growth[str(medium)][len(exponential_growth[str(medium)]) - 1]]]
      time_exp = np.diff(time_exp)
      mu = float(((DO2 - DO1)*2.303) / time_exp)
      Tg = np.log(2)/mu
      results[str(medium)] = [mu, Tg]
      
   return results

def exponential_growth_model(t, No, mu):
   return No * np.exp(mu * t)

## Obtaining Kinetic parameters base on Modelling approach.
def adjusting_exponential_curve(Media, time, exponential_growth, p0):
   
   results = {}
   
   for medium in range(len(Media)):
      t_exp = time[exponential_growth[str(medium)][0]:(exponential_growth[str(medium)][0]+len(exponential_growth[str(medium)]))]
      OD_exp = Media[str(medium)][exponential_growth[str(medium)][0]:(exponential_growth[str(medium)][0]+len(exponential_growth[str(medium)]))]
      params, cv = curve_fit(exponential_growth_model, t_exp, OD_exp , p0)
      No, mu = params
      Tg = np.log(2)/mu
      # determine quality of the fit
      squaredDiffs = np.square(OD_exp - exponential_growth_model(t_exp, No, mu))
      squaredDiffsFromMean = np.square(OD_exp - np.mean(OD_exp))
      rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
      results[str(medium)] = [No, mu, Tg, rSquared]
      
   return results

## Plotting growth curves for all media with exponential adjust.
def plottingCurves(time, Media, scale, time_unit, exponential_growth, parameters_modal):
   
   for medium in range(len(Media)):
      
      plt.plot(time, Media[str(medium)], linestyle='-', marker='o', label = f"Meio {medium}")
      t_exp = time[exponential_growth[str(medium)][0]:(exponential_growth[str(medium)][0]+len(exponential_growth[str(medium)]))]
      OD_exp = Media[str(medium)][exponential_growth[str(medium)][0]:(exponential_growth[str(medium)][0]+len(exponential_growth[str(medium)]))]
      r_square = str(parameters_modal[str(medium)][3])
      r_square = r_square[0:5]
      plt.plot(t_exp, exponential_growth_model(t_exp, parameters_modal[str(medium)][0], parameters_modal[str(medium)][1]), '--', color='black', label=f'Fase exponencial ajustada {r_square}')
      
   plt.xlabel(f'Tempo ({time_unit})')
   plt.ylabel(f'OD600 ({scale})')
   plt.axhline(y=0.5, color='red', linestyle='--', label='Indução com IPT')
   plt.legend()
   plt.yscale(scale)
   plt.savefig(f'growthCurve_{scale}')
   plt.show()
   
def creating_report(file_path, kinetic_parameters, unidade):
   report = open(file_path, 'w')
   for medium in range(len(kinetic_parameters["Matematical"])):
      report.write(f'### Meio {medium} ###\n-Mu em 1/{unidade} (matematico): {str(kinetic_parameters["Matematical"][str(medium)][0])[0:5]}, (modelo exponencial): {str(kinetic_parameters["Exponential_Model"][str(medium)][1])[0:5]};\n-Tg em {unidade} (matematico): {str(kinetic_parameters["Matematical"][str(medium)][1])[0:4]}, (modelo exponencial): {str(kinetic_parameters["Exponential_Model"][str(medium)][2])[0:4]};\n-Ajuste R_quadrado: {str(kinetic_parameters["Exponential_Model"][str(medium)][3])[0:5]}\n\n')
      
def main():
   file_path_input = 'Data.txt' ## it must be in scv formmat, with time in first column and media DO in the others.
   file_path_output = 'Report.txt'
   time_interval = 30
   time_unit = 'min'
   graph_scale = ['linear', 'log']
   faixa = [0.1, 0.2, 0.1, 0.2, 0.25, 0.1, 0.1] ## It must be in the same order of media in Data.tx
   p0 = (0.007, 0.03) # start with values near those we expect to best adjust on exponential phase (No, mu) according with litterature.
   time, Media = getting_object(file_path_input)
   exponential_growth_indexes = finding_exponential_growth(time_interval, Media, time, faixa)
   kinetic_parameters = {}
   kinetic_parameters['Matematical'] = (obtaining_kinetic_parameters_Matematical(Media, time, exponential_growth_indexes))
   kinetic_parameters['Exponential_Model'] = adjusting_exponential_curve(Media, time, exponential_growth_indexes, p0)
   plottingCurves(time, Media, graph_scale[1], time_unit, exponential_growth_indexes, kinetic_parameters['Exponential_Model'])
   creating_report(file_path_output, kinetic_parameters, time_unit)
   
main()

## Make sure to install libs using pip instal and run 'python3 Growth_curve_analysis.py' in terminal to generate files report.txt with all the kinetics data calculated and fitted, and the png image of final growth curve with all media.

## Once the code is run, the files are automatically replaced by the new analysis, the names must be changed for keeping the file or it can be copied to another file.