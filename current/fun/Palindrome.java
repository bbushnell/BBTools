package fun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Advanced string processing utility for finding the longest palindromic substring
 * with configurable mismatch tolerance and optional reverse complement support.
 * Supports processing sequences from files or command-line arguments with
 * biological sequence compatibility through reverse complement matching.
 *
 * @author Brian Bushnell
 */
public class Palindrome {
	
	/**
	 * Main entry point that processes command-line arguments and finds longest
	 * palindrome across all input sequences. Supports file input, direct sequence
	 * input, and configuration parameters for mismatch tolerance and loop size.
	 * @param args Command-line arguments: sequences, filenames, "rcomp"/"rc" flag,
	 * numeric mismatch count, or loop=N parameter
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
	 * Reads sequences from a FASTA-formatted file, parsing sequence headers
	 * and concatenating sequence lines into complete sequences.
	 * @param fname Path to input file containing sequences
	 * @return List of sequences extracted from the file, empty list on error
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
	 * Finds the longest palindromic substring using both even and odd-length
	 * center expansion at each position. Compares lengths from both methods
	 * and returns the longer palindrome found.
	 * @param s Input sequence to analyze for palindromes
	 * @return Longest palindromic substring found in the sequence
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
	 * Finds the longest palindrome allowing for internal loops or bulges up to
	 * maxloop size. Tests both even and odd-length palindromes with varying
	 * loop sizes from 0 to maxloop.
	 * @param s Input sequence to analyze
	 * @param maxloop Maximum internal loop size allowed in palindromes
	 * @return Longest palindrome found with allowed loop tolerance
	 */
	public static String longestPalindrome(String s, int maxloop){
		int longestLength=0;
		int longestStart=0;
		String p="";
		for(int loop=0; loop<=maxloop; loop++) {
			for(int i=0; i<s.length(); i++){
				if((loop&1)==1) {//odd
					int lenOdd=palindromeLengthOdd(s, i, i+loop+1);
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
	 * Calculates length of odd-centered palindrome starting from middle position.
	 * Delegates to the two-parameter version with identical start and end points.
	 * @param s Input sequence
	 * @param middle Center position for palindrome expansion
	 * @return Length of the longest odd-centered palindrome at this position
	 */
	public static int palindromeLengthOdd(String s, int middle){
		return palindromeLengthOdd(s, middle, middle);
	}
	/**
	 * Calculates maximum odd-centered palindrome length between positions a and b,
	 * allowing up to maxMismatches mismatches. Expands outward while tracking
	 * mismatches and stores final boundaries in static variables a_ and b_.
	 * @param s Input sequence
	 * @param a Left starting position
	 * @param b Right starting position
	 * @return Length of palindrome found, 0 if invalid or too short
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
	 * Calculates length of even-centered palindrome starting from middle position.
	 * Delegates to the two-parameter version with adjacent start positions.
	 * @param s Input sequence
	 * @param middle Left center position for palindrome expansion
	 * @return Length of the longest even-centered palindrome at this position
	 */
	public static int palindromeLengthEven(String s, int middle){
		return palindromeLengthEven(s, middle, middle+1);
	}
	/**
	 * Calculates maximum even-centered palindrome length between positions a and b,
	 * allowing up to maxMismatches mismatches. Expands outward while tracking
	 * mismatches and stores final boundaries in static variables a_ and b_.
	 * @param s Input sequence
	 * @param a Left starting position
	 * @param b Right starting position
	 * @return Length of palindrome found
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
	 * Determines if two characters match according to current matching rules.
	 * Uses direct equality when rcomp=false, or reverse complement matching
	 * when rcomp=true for biological sequence analysis.
	 * @param a First character to compare
	 * @param b Second character to compare
	 * @return true if characters match under current rules
	 */
	static boolean matches(char a, char b) {
		return a==(rcomp ? baseToComp[b] : b);
	}
	
	/**
	 * Creates lookup table mapping nucleotide bases to their complements.
	 * Maps A/a to T, T/t/U/u to A, C/c to G, G/g to C, with all other
	 * characters defaulting to N for unknown bases.
	 * @return Character array mapping bases to their complements
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
	
	/** Whether to use reverse complement matching for biological sequences */
	static boolean rcomp=false;
	/** Maximum number of mismatches allowed in palindrome detection */
	static int maxMismatches=0;
	/** Maximum internal loop size allowed in palindromic structures */
	static int maxLoop=0;
	/** Lookup table mapping nucleotide bases to their complements */
	static final char[] baseToComp=makeBaseToComp();
	
	static int a_, b_;
	
}
