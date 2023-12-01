# AuRaCh - Automated Radiation Chemistry Simulations

AuRaCh is a tool for simulation of radiolysis of liquids under exposure to ionizing radiation.
Albeit it has been primarily written for liquid-phase transmission electron microscopy, it can be used for other ionizing sources than electrons, such as x-rays, as well.
Its main advantage is the automation of the creation of coupled differential equation sets.
Consequently, AuRaCh can be employed for non-irradiated  chemical reaction kinetic engineering, as well, simply by setting the dose rate to 0. 

Its working principle is described in:

B. Fritsch, T.S. Zech, M.P. Bruns, A. Körner, S. Khadivianazar, M. Wu, N. Zargar Talebi, S. Virtanen, T. Unruh, M.P.M. Jank, E. Spiecker, A. Hutzler, Radiolysis‐Driven Evolution of Gold Nanostructures –
Model Verification by Scale Bridging _in situ_ Liquid‐Phase Transmission Electron
Microscopy and X‐Ray Diffraction, _Advanced Science_ 2022, 9 (25), 2202803. [DOI:10.1002/advs.202202803](https://doi.org/10.1002/advs.202202803).

If you are new to running python scripts, I propose to use Anaconda: https://www.anaconda.com/products/distribution.
I used Notepad++ to set up the reaction text file, so the spacing is adjusted to display properly using this software. https://notepad-plus-plus.org/ 


Naturally, this code is still under ongoing development. I invite you to contribute :)

## Folder description

The folder "working example" provides code to perform a simulation of the HAuCl<sub>4</sub> reaction set taken off [the original AuRaCh publication](https://doi.org/10.1002/advs.202202803) as exemplary usecase scenario.
To run the code, please acknowledge the required dependencies.
Its default environment should be able to run all files uploaded here.


"AuRaCh reaction networks - Radiolytic Acidity" contains the reaction networks and initial concentrations/g-values file for the simulations performed in 

B. Fritsch, A. Körner, T. Couasnon, R. Blukis, M. Taherkhani, L.G. Benning, M.P.M. Jank, E. Spiecker, A. Hutzler, Tailoring the Acidity of Liquid Media with Ionizing Radiation: Rethinking the Acid-Base Correlation Beyond pH, _The Journal of Physical Chemistry Letters_, 2023, 14 (20), 4644–4651. [DOI:10.1021/acs.jpclett.3c00593](https://doi.org/10.1021/acs.jpclett.3c00593).

A reaction file for simulating an aqueous solution containing Fe(II)/Fe(III) can be found in "AuRaCh reaction network - Fe". It is introduced in

T. Couasnon, B. Fritsch, M.P.M. Jank, R. Blukis, A. Hutzler, L. Benning, Goethite Mineral Dissolution to Probe the Chemistry of Radiolytic Water in Liquid-Phase Transmission Electron Microscopy, _Advanced Science_ 2023, 10 (25), 2301904. [DOI:10.1002/advs.202301904](https://doi.org/10.1002/advs.202301904).

To create subsets of reaction networks (sparse reaction sets), the script in the folder "Sparse Reaction Files" can be used. Alongside, you'll find three examples. Those, as well as the concept, are explained in 

B. Fritsch, P. Malgareti, J. Harting, Karl J.J. Mayrhofer, A. Hutzler, Precision of Radiation Chemistry Networks: Playing Jenga with Kinetic Models for Liquid-Phase Electron Microscopy,
_Precision Chemistry_, [DOI:10.1021/prechem.3c00078](https://doi.org/10.1021/prechem.3c00078).


Please use the AuRaCh code file or the GUI to run the simulations or contact me in case of further questions.
