# MRI_IR_Calculator_MatLab

## MRI (D)IR TI/T1++ calculator for MatLab, by Øystein Bech-Aase (formerly Bech Gadmar)
### Based on original code by Atle Bjørnerud, in collaboration with Wibeke Nordhøy

Just run the main .m file in Matlab. Tested mostly on R2018-R2019.
  
TrueT2-DIR is based on Madhuranthakam et al, Mag Res Med 2012 67(1):81-88.
  
For more info, see comments in the IR_Calculator.m file.

## DIR T1 vs TI calculator for LibreOffice Calc by GadOE

With this spreadsheet tool you can calculate TrueT2-DIR TI sequence parameters 
based on measured T1 values, for a predetermined set of other parameters.
  
The T1-vs-TI IR-Calc tool could be deployed in clinical settings without MatLab access. 
Simulations and regressions therein have to be redone for each clinical setting.

## Online Resources

- Miho Kita's online [tool for calculating TI for various sequence types](https://seichokai.jp/fuchu/null_point_english/)
- XRayPhysics' explanation of [simple inversion mathematics](https://xrayphysics.com/contrast.html)
- On T1 times in tissues: [Lalande et al, MRI 2016](https://www.sciencedirect.com/science/article/pii/S0730725X16301266)


---
Have fun inverting!
_Øystein_


## Addendum: Using IR-Calc for non-brain applications

The user T1/T2 times can be set to any tissues you like, as long as you figure out what acts as "WM/GM/MS/CSF" for your purposes.

These are default T1/T2 times for the "1.5 T" set of user relaxation times, taken from GE's values and literature:
- ?? :  T1      T2      (ms)
- WM :  660     80
- GM : 1200     95
- MS : 1200    100
- CSF: 4308   3500

Here are some experimental values used for imaging synovitis in the knee in one of our studies:
- WM : 1000     50      for healthy synovium
- GM : (not used)
- MS : 1100     80      for synovitis
- CSF: 3600   1200      for synovial fluid

These values are mostly guesswork for now, in part based on attempted measurement of said values.
But the point is, IR-Calc can be used this way. Of course, it can be a challenge to procure good value estimates.
