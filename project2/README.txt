AUTHOR: Paul Villanueva

USE:
Put bcb567_project2_classes, bcb567_project2_utils, and bcb567_project2_main in the same directory.  Then call bcb567_project2_main by typing 

	python bcb567_project1_main file1 file2 wlcut
	
into the command line, where file1 is a .txt file in FASTA format, file2 is a .txt file containing a single line of 1s and 0s, wlcut is a positive integer.

The output is a printout to the console in the following format:


                              Word model: x
                                   wlcut: x 
    Number of positions in largest block: x
            Positions in superword array: x  x  ... x
              Positions in largest block: x  x  ... x
              Superword of largest block: x
	
where each of the xs are stand-ins for the appropriate values.  

As an example, entering "python bcb567_project1_main.py D.short.txt word_model.txt 3" will produce the output associated with the sample superword array given on blackboard.

FILES:
All files written and compiled on Python 2.7.14.

bcb567_project1_classes: contains the SuperwordArray class.
	LocalAlignment - takes in a dna sequence, a word model, and a wlcut and computes the superword array and the largest block.
	
bcb567_project1_utils: contains the methods called in bcb567_project1_classes.    See documentation in each file for further 
information.

bcb567_project1_main: contains the driver for the entire package.  When called as indicated in the USE section, the driver first 
reads the two .txt files to create a dna string and a word model, then creates a SuperwordArray object initialized with these two strongs and wlcut.  It then prints out the information outlined above.

bcb567_project2_sample_output is a screen shot of running "python bcb567_project1_main.py D.short.txt word_model.txt 3" on the command line.

