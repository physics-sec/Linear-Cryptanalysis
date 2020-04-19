package main

import "fmt"

var NUM_P_C_PAIRS int = 5000
var SBOX_BITS int = 4
var NUM_SBOXES int = 4
var NUM_ROUNDS int = 4
var MIN_BIAS float32 = 0.008
var MAX_BLOCKS_TO_BF int = 3

func apply_mask(value, mask int)int {
    // retrieve the parity of mask/value
    interValue := value & mask
    total := 0
    for interValue > 0 {
        temp := interValue & 1
        interValue = interValue >> 1
        if temp == 1 {
            total = total ^ 1
        }
    }
    return total
}

// convert a list of bits to an integer
func bits_to_num(inputbits []int) int {
    Y_input := 0
    for _, input_pos := range inputbits {
        Y_input = Y_input | (1 << (SBOX_BITS - input_pos))
    }
    return Y_input
}

// convert an integer to a list of it's bits
func num_to_bits(num int) []int {
    var bits []int
    for index := 0; index < SBOX_BITS; index++ {
        if (1 << index) & num > 0 {
            bits = append(bits, SBOX_BITS - index)
        }
    }
    return bits
}

func do_sbox(num int) int {
    s := map[int]int{0:0xE, 1:0x4, 2:0xD, 3:0x1, 4:0x2, 5:0xF, 6:0xB, 7:0x8, 8:0x3, 9:0xA, 0xA:0x6, 0xB:0xC, 0xC:0x5, 0xD:0x9, 0xE:0x0, 0xF:0x7}
    return s[num]
}

/*

func create_bias_table() []float64 {
    tablesize := 1 << SBOX_BITS
    var bits []int

    var table [][3]float64

    for x := 1; x < tablesize; x++ {
        for y := 1; y < tablesize; y++ {
            matches := 0
            for num := 0; num < tablesize; num++ {
                in_mask  := apply_mask(num, x)
                out_mask := apply_mask(do_sbox(num), y)
                if in_mask == out_mask {
                    matches++
                }
            }
       
            bias := (matches / tablesize) - 1/2
            if bias > 0 {
                table = append(table, [3]float64{x, y, bias} )
            }
        }
    }
 }
*/



/*
def create_bias_table():
    if not initialized: exit('initialize the library first!')

    tablesize = 1 << SBOX_BITS
   
    table = []
    for x in range(1, tablesize):
        for y in range(1, tablesize):
            matches = 0
            for num in range(tablesize):
                # calculate the parity of the number before going in the sbox
                in_mask  = apply_mask(num, x)
                # calculate the parity of the number after going out the sbox
                # do_sbox is supposed to perform the substitution, make sure is well defined!
                out_mask = apply_mask(do_sbox(num), y)
                # if the parity is the same in both cases, add a match
                if in_mask == out_mask:
                    matches += 1

            # calculate the bias
            bias = (matches / tablesize) - 1/2
            # if the bias is greater than 0, save the 'x' and 'y' combination
            if bias > 0:
                table.append( [x, y, bias] )

    # return the table of biases that are greater than 0
    return table
*/

func main() {
    fmt.Println(do_sbox(1))
    create_bias_table()
    fmt.Println("jejej")
}
