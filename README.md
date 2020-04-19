# Linear-Cryptanalysis


## What is this?

This repository presents a linear cryptanalysis library that tries to break SPN ciphers in a fully automatic way.  

Right now it only supports SPNs with just one type of sbox, but extending it to support multiple types of sbox should be easy.


## How does it work?

Glad you asked, the algorithm is very simple.  

### First

We compute the bias table. There is nothing fancy about this.  
All the biases over 0 are kept and sorted in descending order.

### Second

In an SPN, we have C sboxes in R rows (one per round).  

For each sbox in the first row, we calculate all possible 'moves' (based on the bias table).  
This will give us *C\*len(bias table)* possible first steps, each asociated with a bias.  
Each step will reach different 'input bits' in some of the sboxs in the second row.  

For each possible first step, we do the following:  
We now position ourselves at the sboxs of the second row, the ones that were reached by our first sbox.  
For each one of these, we calculate all the possible steps, keeping in mind that we must choose the biases with the same 'input bits'.  
For example, if the sbox of the first step reached into the second and fourth bit of an sbox in the second layer, then we can only choose biases that have 2 and 4 as input for that particular sbox.  

Once we calculated all the possible moves for all the sboxs in the second layer, we combine all this moves in all possible ways.  
For example, if in the second layer, 3 sboxs were reached, and each can move in 4 different ways, there are 64 (4x4x4) combinations.  
For each possible combination, we do the same as in the first step, we calculate which sboxes in the third layer were reached and in which bits.  
Now we position ourselves at the sboxes in the third row and continue the process until the last row.  

At the end, we return all possible paths (or linear approximations) to the last row, each associated with a bias, calculated according to the Piling-Up Lemma.  


### Third

We choose one linear approximation from the previous step (normally the one with the highest bias).  

Now, we generate multiple plaintext/ciphertext pairs and use the chosen linear approximation to obtain bits of the last round key. (again, nothing fancy here, just normal linear cryptanalysis.)

## How to use

### Example  
There is an example in the `break-basic_SPN.py` file.  

### Functions you may want to use

`initialize()`  
Initializes the library with the SPN's properties.  

`create_bias_table()`  
Creates the bias table for the sbox.  

`get_linear_aproximations(bias_table)`  
Returns all possible linear approximations (that respect the MAX_BLOCKS_TO_BF filter).  

`analize_cipher()`  
Creates the bias table (calling `create_bias_table`) and sorts it (if it is longer than 1000 rows, keeps the best 1000).  
Then, calls `get_linear_aproximations` and sorts the results (deleting the approximations that have a bias lower than MIN_BIAS) and returns the sorted list of approximations.  

`get_biases(_c_pairs, linear_aproximation)`  
Returns a list of of biases. The index of the bias is the key used to obtain such bias.

## Considerations

Keep in mind that you might use multiple linear approximations to recover different bits of the last round key.  

If the bias table is bigger than 1000 rows, the library keeps the best 1000, you might want to proceed differently.  

## More information

To learn about linear cryptanalysis, read [this](https://www.engr.mun.ca/~howard/PAPERS/ldc_tutorial.pdf) awesome paper by Howard M. Heys and read *Modern Cryptanalysis: Techniques for Advanced Code Breaking* by Christopher Swenson.


## Credit

Thanks to [hkscy](https://github.com/hkscy/Basic-SPN-cryptanalysis) for the great *Basic SPN* implementation.
