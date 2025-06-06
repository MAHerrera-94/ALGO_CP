# ALGO_CP
A bespoke sequence-generating algorithm for the design of novel acyl carrier proteins.
The following code allows the user to generate new ACP sequences from a high-quality multiple sequence alignment of homologoues. Here we have provided the input AcpP MSA used in our studies. For ease of execution, ensure the file "AcpP_input_MSA.fasta" file is in the same directory as ALGO_CP.py. If desired, this can be specified explicitly in line 14.
ALGO-CP uese two complementary sub-algorithms - route-AC and route AP - to create novel ACP-like sequences. The balance of these two sub-algorithms can be configured by a weighting coefficient r (line 120).
The physicochemical weights used by route-AP can be changed in lines 112-114.
