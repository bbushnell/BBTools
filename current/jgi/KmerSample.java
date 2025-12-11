package jgi;

import fileIO.TextFile;

/**
 * K-mer sampling and storage using bit arrays for memory-efficient representation.
 * Reads sequences from FASTA files and stores k-mer presence information in compact bit arrays.
 * Each k-mer is represented by a single bit, allowing storage of large k-mer sets
 * in minimal memory.
 *
 * @author Brian Bushnell
 * @date Oct 10, 2012
 */
public class KmerSample {
	
	
	/**
	 * Builds a bit array marking presence of all k-mers (length K) found in a FASTA file.
	 * Each k-mer maps to a single bit in the returned array; ambiguous bases reset the rolling k-mer.
	 * @param K K-mer length (<=31)
	 * @param filename FASTA file path
	 * @return Bit array marking observed k-mers
	 */
	public static int[] makeKmerSet(int K, String filename){
		
		//Number of bits in a kmer
		int kbits=2*K;
		
		//Number of possible kmers
		long kmerSpace=(1L<<kbits);
		
		//Make an array of the correct size, remembering that int is 32 bits
		int[] array=new int[(int)(kmerSpace/32)];
		
		//Current kmer
		long kmer=0;
		
		//Length of current kmer
		int len=0;
		
		//This will create a bitmask of 00000...0000111111...11111, where the number if 1's is equal to kbits.
		long mask=(kbits>63 ? -1L : ~((-1L)<<kbits));
		
		//Initialize an input stream for the fasta file
		TextFile tf=new TextFile(filename, false);
		
		//Grab the first line of the fasta file
		String line=tf.nextLine();
		
		while(line!=null){
			
			if(line.length()<1){
				//The line is empty, so ignore it (should never happen in a proper fasta file)
			}else if(line.charAt(0)=='>'){
				//The line is name of a new contig/scaffold, so reset the kmer
				kmer=0;
				len=0;
			}else{
				//Otherwise, generate kmers
				
				for(int i=0; i<line.length(); i++){
					
					//The base at location "i" in the string
					char letter=line.charAt(i);
					
					//The 2-bit numeric code for the base
					int code;
					
					if(letter=='A'){code=0;}
					else if(letter=='C'){code=1;}
					else if(letter=='G'){code=2;}
					else if(letter=='T'){code=3;}
					else{code=-1;}
					
					if(code<0){
						//The base was an N or degenerate letter, so reset the kmer
						kmer=0;
						len=0;
					}else{
						//insert the code into the current kmer
						kmer=(kmer<<2); //left shift by 2
						kmer=(kmer|code); //or with the code
						kmer=(kmer&mask); //and with the mask to prevent going past the intended kmer length
						len++; //Increment the length of the kmer
						
						if(len>=K){
							//If the kmer is long enough, then add it to the array
							
							//The index in the array is the upper bits of the kmer.  Each location in the array is 32 bits.
							int index=(int)(kmer/32);
							
							//The bit within the word of the array is the lower 5 bits of the kmer
							int bit=(int)(kmer%32);
							
							//A bitmask to set the correct bit in the array to 1.
							int x=(1<<bit);
							
							//OR the array location with the new mask.
							array[index]=(array[index] | x);
						}
					}
				}
			}
			
			//Grab the next line
			line=tf.nextLine();
		}
		
		//Close your input stream
		tf.close();
		
		return array;
	}
	
	/**
	 * Tests whether a specific k-mer (encoded as a long) is present in the bit array produced by makeKmerSet.
	 * @param kmer Encoded k-mer
	 * @param array Bit array from makeKmerSet
	 * @return true if the k-mer bit is set
	 */
	public static boolean containsKmer(long kmer, int[] array){
		
		//The index in the array is the upper bits of the kmer.  Each location in the array is 32 bits.
		int index=(int)(kmer/32);
		
		//The bit within the word of the array is the lower 5 bits of the kmer
		int bit=(int)(kmer%32);
		
		//A bitmask to test the correct bit in the array to 1.
		int x=(1<<bit);
		
		if((array[index]&x)==0){//Check to see if the bit is set in the array
			return false;
		}else{
			return true;
		}
		
	}
	
}
