package fun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Advanced palindrome detection utility for identifying the longest palindromic
 * substring with configurable mismatch tolerance, loop support, and optional
 * reverse complement matching for biological sequences.
 * Supports both file input and direct string processing.
 *
 * @author Brian Bushnell
 */
public class Palindrome2 {
	
	/**
	 * Program entry point for palindrome detection.
	 * Processes command-line arguments for configuration and input sequences.
	 * Arguments can include: rcomp/rc (reverse complement mode), numeric values
	 * (max mismatches), loop=N (max loop size), file paths, or direct sequences.
	 *
	 * @param args Command-line arguments for configuration and input
	 */
	public static void main(String[] args){
		
		ArrayList<String> sequences=new ArrayList<String>();
		
		for(String s : args) {
			if(s.equalsIgnoreCase("rcomp") || s.equalsIgnoreCase("rc")) {
				rcomp=true;
			}else if(Character.isDigit(s.charAt(0))) {
				maxMismatches=Integer.parseInt(s);
			}else if(s.startsWith("loop")) {
				maxLoop=Integer.parseInt(s.split("=")[1]);
			}else if(new File(args[0]).exists()) {
				sequences.addAll(getSequence(args[0]));
			}else {
				sequences.add(s);
			}
		}
		
		String longest="";
		for(String s : sequences) {
			String p;
			if(maxLoop<1) {
				p=longestPalindrome(s);
			}else {
				p=longestPalindrome(s, maxLoop);
			}
			if(p.length()>longest.length()) {
				longest=p;
			}
		}
		System.out.println("Longest palindrome is length "+longest.length()+":\n'"+longest+"'");
	}
	
	/**
	 * Reads sequences from a FASTA-like file format.
	 * Parses lines starting with '>' as headers and concatenates subsequent
	 * lines as sequence data until the next header is encountered.
	 *
	 * @param fname Path to the input file containing sequences
	 * @return List of sequence strings extracted from the file
	 */
	public static ArrayList<String> getSequence(String fname){
		ArrayList<String> list=new ArrayList<String>();
		try {
			final BufferedReader reader=new BufferedReader(new FileReader(fname));

			StringBuilder sb=new StringBuilder();
			for(String line=reader.readLine(); line!=null; line=reader.readLine()) {
				if(line.length()>0) {
					if(line.charAt(0)=='>'){
						if(sb.length()>0) {
							list.add(sb.toString());
							sb.setLength(0);
						}
					}else{
						sb.append(line);
					}
				}
			}
			if(sb.length()>0) {list.add(sb.toString());}
			reader.close();
		}catch(Exception e){
			
		}
		return list;
	}
	
	/**
	 * Finds the longest palindromic substring in the input string.
	 * Checks both odd-length (single center) and even-length (dual center)
	 * palindromes at each position. Uses mismatch tolerance if configured.
	 *
	 * @param s Input string to analyze for palindromes
	 * @return The longest palindromic substring found
	 */
	public static String longestPalindrome(String s){
		int longestLength=0;
		int longestStart=0;
		String p="";
		for(int i=0; i<s.length(); i++){
			int lenEven=palindromeLengthEven(s, i);
			if(lenEven>longestLength){
				longestLength=lenEven;
				longestStart=i-lenEven/2+1;
				int a=longestStart, b=longestStart+longestLength;
				p=s.substring(a, b);
			}
			int lenOdd=palindromeLengthOdd(s, i);
			if(lenOdd>longestLength){
				longestLength=lenOdd;
				longestStart=i-lenOdd/2;
				p=s.substring(longestStart, longestStart+longestLength);
			}
		}
		return p;
	}
	
	/**
	 * Finds the longest palindromic substring with configurable loop support.
	 * Allows for internal loops/gaps within palindromes by testing different
	 * loop sizes from 0 to maxloop. Alternates between odd and even palindrome
	 * detection based on loop parity.
	 *
	 * @param s Input string to analyze for palindromes
	 * @param maxloop Maximum internal loop/gap size allowed
	 * @return The longest palindromic substring found with loop support
	 */
	public static String longestPalindrome(String s, int maxloop){
		int longestLength=0;
		int longestStart=0;
		String p="";
		for(int loop=0; loop<=maxloop; loop++) {
			for(int i=0; i<s.length(); i++){
				if((loop&1)==1) {//odd
					int lenOdd=palindromeLengthOdd(s, i, i+loop);
					if(lenOdd>longestLength){
						longestLength=lenOdd;
//						longestStart=i-lenOdd/2;
//						p=s.substring(longestStart, longestStart+longestLength-loop);
						p=s.substring(a_, b_+1);
					}
				}else {//even
					int lenEven=palindromeLengthEven(s, i, i+1+loop);
					if(lenEven>longestLength){
						longestLength=lenEven;
//						longestStart=i-lenEven/2+1;
//						int a=longestStart, b=longestStart+longestLength-loop;
//						p=s.substring(a, b);
						p=s.substring(a_, b_+1);
					}
				}
			}
		}
		return p;
	}
	
	/**
	 * Calculates the length of an odd-length palindrome centered at the
	 * specified position. Delegates to the two-parameter version with
	 * identical start and end positions.
	 *
	 * @param s Input string to analyze
	 * @param middle Center position of the potential palindrome
	 * @return Length of the palindrome found, or 0 if none
	 */
	public static int palindromeLengthOdd(String s, int middle){
		return palindromeLengthOdd(s, middle, middle);
	}
	/**
	 * Calculates the length of an odd-length palindrome between positions a and b.
	 * Expands outward from the initial positions, allowing for mismatches up to
	 * the configured threshold. Updates global position variables a_ and b_.
	 *
	 * @param s Input string to analyze
	 * @param a Starting left position
	 * @param b Starting right position
	 * @return Length of the palindrome found, accounting for mismatches
	 */
	public static int palindromeLengthOdd(String s, int a, int b){
		int length=b-a-1;
		int mismatches=0;
		while(a>=0 && b<s.length() && mismatches<=maxMismatches){
			if(!matches(s.charAt(a), s.charAt(b))){
				mismatches++;
				if(mismatches>maxMismatches) {break;}
			}
			length+=2;
			a--;
			b++;
		}
		if(a<0 || b>=s.length() || mismatches>maxMismatches) {a++; b--;}
		a_=a;
		b_=b;
		if(length==1 && rcomp && maxMismatches<1) {return 0;}
		return length<0 ? 0 : length;
//		return b-a+1;
	}

	/**
	 * Calculates the length of an even-length palindrome centered between
	 * the specified position and the next position. Delegates to the
	 * two-parameter version with adjacent start positions.
	 *
	 * @param s Input string to analyze
	 * @param middle Left center position of the potential palindrome
	 * @return Length of the palindrome found, or 0 if none
	 */
	public static int palindromeLengthEven(String s, int middle){
		return palindromeLengthEven(s, middle, middle+1);
	}
	/**
	 * Calculates the length of an even-length palindrome between positions a and b.
	 * Expands outward from the initial positions, allowing for mismatches up to
	 * the configured threshold. Updates global position variables a_ and b_.
	 *
	 * @param s Input string to analyze
	 * @param a Starting left position
	 * @param b Starting right position
	 * @return Length of the palindrome found, accounting for mismatches
	 */
	public static int palindromeLengthEven(String s, int a, int b){
		int length=b-a-1;
		int mismatches=0;
		while(a>=0 && b<s.length() && mismatches<=maxMismatches){
			if(!matches(s.charAt(a), s.charAt(b))){
				mismatches++;
				if(mismatches>maxMismatches) {break;}
			}
			length+=2;
			a--;
			b++;
		}
		if(a<0 || b>=s.length() || mismatches>maxMismatches) {a++; b--;}
		a_=a;
		b_=b;
		return length;
//		return b-a+1;
	}
	
	/**
	 * Determines if two characters match based on current matching mode.
	 * In reverse complement mode, compares character a with the complement
	 * of character b using the base complement lookup table.
	 *
	 * @param a First character to compare
	 * @param b Second character to compare
	 * @return true if characters match according to current mode
	 */
	static boolean matches(char a, char b) {
		return a==(rcomp ? baseToComp[b] : b);
	}
	
	/**
	 * Creates a lookup table for DNA base complements.
	 * Maps A/a to T, C/c to G, G/g to C, T/t/U/u to A, with all other
	 * characters mapped to 'N'. Used for reverse complement matching.
	 * @return Character array mapping ASCII values to their complements
	 */
	static char[] makeBaseToComp() {
		char[] array=new char[128];
		Arrays.fill(array, 'N');
		array['A']=array['a']='T';
		array['C']=array['c']='G';
		array['G']=array['g']='C';
		array['T']=array['t']=array['U']=array['u']='A';
		return array;
	}
	
	private static class Drome{
		
		/**
		 * Constructs a Drome object representing a palindromic sequence.
		 * Extracts the substring from the reference and begins analysis
		 * of the palindrome structure by examining character pairs.
		 *
		 * @param ref Reference string containing the palindrome
		 * @param start Starting position in the reference string
		 * @param stop Ending position in the reference string
		 */
		public Drome(String ref, int start, int stop) {
			s=ref.substring(start, stop+1);
			final int half=s.length()/2;
			for(int i=0, j=s.length()-1; i<half; i++, j--) {
				final char ci=s.charAt(i);
				final char cj=s.charAt(j);
			}
		}
		
		/** The palindromic sequence string */
		String s;
		/** Size of internal loop/gap within the palindrome */
		int loop;
		/** Number of mismatched positions in the palindrome */
		int mismatches;
	}
	
	/** Whether to use reverse complement matching for biological sequences */
	static boolean rcomp=false;
	/** Maximum number of mismatches allowed in palindrome detection */
	static int maxMismatches=0;
	/** Maximum internal loop/gap size allowed in palindromes */
	static int maxLoop=0;
	/** Lookup table mapping DNA bases to their complements */
	static final char[] baseToComp=makeBaseToComp();
	
	/** Global variable storing the right boundary of the last found palindrome */
	/** Global variable storing the left boundary of the last found palindrome */
	static int a_, b_;
	
}
