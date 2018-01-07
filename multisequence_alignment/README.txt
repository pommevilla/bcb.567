AUTHOR: Paul Villanueva

USE:
Put bcb567_project3_classes, bcb567_project3_utils, and bcb567_project3_main in the same directory.  Then call bcb567_project3_main by typing 

	python bcb567_project3_main seqs_file word_model wlcut
	
into the command line, where seqs_file is a .txt file in FASTA format, word_model is a .txt file containing a single line of 1s and 0s, and wlcut is a positive integer.

The output is a printout to the console in the following format:

>Seq1
xxxx...
>Seq2
NANACNCGNGNNNTNTTNT
xxx...
>Seqn

Word model: xxx...
wlcut: x

The length of a longest chain of superword blocks: xx...
The number of superword blocks in the chain: x...

Block 1
Seq1       x        xxx
Seq2       x        xxx
...
Seqn       x        xxx      

...

Block n 
...	Block 1
Seq1       x        xxx
Seq2       x        xxx
...
Seqn       x        xxx

where each of the xs are stand-ins for the appropriate values.  

As an example, entering "python bcb567_project3_main.py multi_2.txt wm_101.txt 2" will produce the output associated with the sample superword array given on blackboard.

FILES:
All files written and compiled on Python 2.7.14.

bcb567_project3_classes: contains the MultiSequenceAlignment and SuperwordArray classes.
	SuperwordArray - Unchanged from assignment #2.
    MultiSequenceAlignment - takes in a list of dna strings, a word model, and a wlcut and reports a longest chain of blocks between the strings.
	
bcb567_project3_utils: contains the methods called in bcb567_project3_classes.    See documentation in each file for further 
information.

bcb567_project3_main: contains the driver for the entire package.  When called as indicated in the USE section, the driver first 
reads the two .txt files to create a list of dna strings and a word model, then creates a MultiSequenceAlignment object initialized with these two strings and wlcut.  It then prints out the information outlined above.

bcb567_project3_sample_output is a screen shot of running "python bcb567_project3_main.py D.short.txt word_model.txt 3" on the command line.

