Installation of MMB Chimera extension
=====================================

Mac OSX
--------
Install Chimera 32 bits. MMB doesn't currently work in 64 bits on Mac OSX, sorry.  
Extract MMB_chimera_extension_OSX.32bits.tgz.  
Copy all dylib files from **libs** to Chimera **lib** directory (usually /Applications/Chimera.app/Contents/Resources/lib)  

Run Chimera.
Go to the "Tools" prefrences windows : Favorites -> Add to Favorites/Toolbar...  
In the bottom panel, click on Add... and choose the directory where you extracted the archive. It should contain MMB_UI and pyMMB directories.  
Click on Save.  

MMB GUI should be available under Tools -> MMB -> MMB User Interface  

Windows
--------
Install Chimera.  
Extract MMB_chimera_extension_Windows.zip  

Run Chimera.  
Go to the "Tools" prefrences windows : Favorites -> Add to Favorites/Toolbar...  
In the bottom panel, click on Add... and choose the directory where you extracted the archive. It should contain MMB_UI and pyMMB directories.  
Click on Save.  

MMB GUI should be available under Tools -> MMB -> MMB User Interface  


Quick Tutorial
=====================================
We are going to generate a small RNA hairpin, as in the first tutorial for MMB in command line.
The main difference is that we don't have stages anymore (a similar feature will be available in the future).

Open Chimera and the MMB user interface.

You should be on the Input tab.

Fill the widgets for a new RNA strand with the following chain ID, first residue number and sequence:
RNA A 2656 UACGUAAGGA
Click on the Add button

We are not going to add more chains, so we can click on the Load button, at the bottom.
Once loaded, you can't change the chains' information anymore.

You should see an RNA strand in the Chimera window, in a beta conformation.

Now, go to the Base Interactions tab.
Here you can specify Base Pair interactions between nucleotides.
Create the following base pairs, either one by one or using the Nucleic Acid Duplex form.
baseInteraction A 2656 WatsonCrick A 2665 WatsonCrick Cis 
baseInteraction A 2657 WatsonCrick A 2664 WatsonCrick Cis 
baseInteraction A 2658 WatsonCrick A 2663 WatsonCrick Cis

You should see dashed lines on the model, indicating the base pairs we want to form during the simulation.

By default everything is flexible and this what we want. So let's skip the Mobilizers tab for now.

However, we need to avoid clashes between atoms during the simulation.
Go to the Contacts tab, and add a contact for all heavy atoms of the chain A.
Don't forget to click on Add !

Let's try to simulate this now.
Go to the simulation tab.
Set the Reporting Intervals to 4.0, the temperature to 30, the dutyCycle to 0.9 and the scrubberPeriod to 40.
Look at the MMB reference guide for an explanation of these parameters.
Check setInitialVelocities to add some randomness at the start of the simulation.
Set the Force Multiplier to 200.

Now, click on Init to initialize the simulation. Warning !!! No parameters can be changed after that !!!
A MD Movie window should appear, we are going to use it later.
Set the number of Reporting Interval (just right near the Init button) at 50 and click on Run.
You should see the structure moving, according to what MMB is simulating.
You should observe that at some point, all the base pairs are satisfied and the structure is not evolving so much. That means that the simulation converged. You can wait until all reporting intervals have been computed or click the Stop button.
Now you can review the trajectory using the MD Movie window.
You can save the resulting structure and the trajectory via Chimera.
If you are not satisfied, you can do another run (you can change the number of reporting intervals between each run). It will start from the last computed structure and the frames will be added to the trajectory.


Now you can try to start over and change the base pairs to generate a GNRA tetraloop from this sequence.
Hints:
#baseInteraction A 2662 Hoogsteen  A 2659 SugarEdge Trans  

#baseInteraction A 2658 Stacking A 2659 Stacking Cis
#baseInteraction A 2660 Stacking A 2661 Stacking Cis
#baseInteraction A 2661 Stacking A 2662 Stacking Cis

If you want a smoother trajectory, you can set the reportingIntervals parameter to 1.0 or less. You will need to increase the number of reporting intervals to reach the same simulation length though.

Unfortunately, the only way to restart an MMB sesssion right now is to quit Chimera and relaunch it.


