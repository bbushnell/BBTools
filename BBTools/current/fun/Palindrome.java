package fun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

public class Palindrome {
	
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
	
	public static int palindromeLengthOdd(String s, int middle){
		return palindromeLengthOdd(s, middle, middle);
	}
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

	public static int palindromeLengthEven(String s, int middle){
		return palindromeLengthEven(s, middle, middle+1);
	}
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
	
	static boolean matches(char a, char b) {
		return a==(rcomp ? baseToComp[b] : b);
	}
	
	static char[] makeBaseToComp() {
		char[] array=new char[128];
		Arrays.fill(array, 'N');
		array['A']=array['a']='T';
		array['C']=array['c']='G';
		array['G']=array['g']='C';
		array['T']=array['t']=array['U']=array['u']='A';
		return array;
	}
	
	static boolean rcomp=false;
	static int maxMismatches=0;
	static int maxLoop=0;
	static final char[] baseToComp=makeBaseToComp();
	
	static int a_, b_;
	
}
