package server;

import java.util.BitSet;
import java.util.HashMap;

import structures.ByteBuilder;

public class PercentEncoding {
	
	public static boolean containsSpecialSymbol(String s){
		if(s==null){return false;}
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			if(isSpecial.get(c)){
//				System.err.print("b");
				return true;}
		}
//		System.err.print("c");
		return false;
	}
	
	public static boolean containsCommonSymbol(String s){
		if(s==null){return false;}
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			if(isCommon.get(c)){return true;}
		}
		return false;
	}
	
	public static String symbolToCode(String s){
//		System.err.print("a");
		if(!containsSpecialSymbol(s)){return s;}
//		System.err.print("d");
		ByteBuilder bb=new ByteBuilder();
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			String code=symbolToCodeArray[c];
			if(code!=null){
//				System.err.print("e("+code+")");
				bb.append(code);
			}else{
//				System.err.print("f");
				bb.append(c);
			}
		}
//		System.err.println("g");
//		System.err.println(bb);
		return bb.toString();
	}
	
	public static String commonSymbolToCode(String s){
		if(!containsCommonSymbol(s)){return s;}
		ByteBuilder bb=new ByteBuilder();
		for(int i=0, max=s.length(); i<max; i++){
			char c=s.charAt(i);
			if(isCommon.get(c)){
				String code=symbolToCodeArray[c];
				assert(code!=null);
				bb.append(code);
			}else{
				bb.append(c);
			}
		}
		return bb.toString();
	}
	
	private static int parseCode(String s, int start){
		if(s==null || start+2>=s.length()){return -1;}
		assert(s.charAt(start)=='%');
		int sum=0;
		for(int i=start+1; i<=start+2; i++){
			sum=sum<<4;
			final char c=s.charAt(i);
			if(c>='0' && c<='9'){
				sum=sum+(c-'0');
			}else if(c>='A' && c<'F'){
				sum=sum+(10+c-'A');
			}else{
				return -1;
			}
		}
		return sum;
	}
	
	public static String codeToSymbol(String s){
		int idx=s.indexOf('%');
		if(idx<0){return s;}
		
		ByteBuilder bb=new ByteBuilder(s.length());
		for(int i=0; i<s.length(); i++){
			char c=s.charAt(i);
			if(c=='%'){
				int sym=parseCode(s, i);
				if(sym<0){bb.append(c);}
				else{
					bb.append((char)sym);
					i+=2;//Skip next 2 characters
				}
			}else{bb.append(c);}
		}
		return (bb.length()==s.length() ? s : bb.toString());
	}
	
	private static HashMap<String, String> makeCodeToSymbolMap() {
		HashMap<String, String> map=new HashMap<String, String>(129);
		assert(reservedSymbol.length==reservedCode.length);
		assert(commonSymbol.length==commonCode.length);
		for(int i=0; i<reservedSymbol.length; i++){
			map.put(reservedCode[i], reservedSymbol[i]);
		}
		for(int i=0; i<commonSymbol.length; i++){
			map.put(commonCode[i], commonSymbol[i]);
		}
		return map;
	}
	
	private static HashMap<String, String> makeSymbolToCodeMap() {
		HashMap<String, String> map=new HashMap<String, String>(257);
		assert(reservedSymbol.length==reservedCode.length);
		assert(commonSymbol.length==commonCode.length);
		for(int i=0; i<reservedSymbol.length; i++){
			map.put(reservedSymbol[i], reservedCode[i]);
		}
		for(int i=0; i<commonSymbol.length; i++){
			map.put(commonSymbol[i], commonCode[i]);
		}
		return map;
	}
	
	private static String[] makeSymbolToCodeArray() {
		final String[] array=new String[128];
		for(int i=0; i<reservedSymbol.length; i++){
			String s=reservedSymbol[i];
			String c=reservedCode[i];
			array[s.charAt(0)]=c;
		}
		for(int i=0; i<commonSymbol.length; i++){
			String s=commonSymbol[i];
			String c=commonCode[i];
			array[s.charAt(0)]=c;
		}
		return array;
	}
	
	private static final BitSet makeBitSet(String[]...matrix){
		BitSet bs=new BitSet(128);
		for(String[] array : matrix){
			for(String s : array){
				char c=s.charAt(0);
				bs.set(c);
			}
		}
		return bs;
	}

	//See https://en.wikipedia.org/wiki/Percent-encoding
	public static final String[] reservedSymbol=new String[] {
		"!", "#", "$", "&", "'", "(", ")", "*", "+", ",", "/", ":", ";", "=", "?", "@", "[", "]"
	};
	
	public static final String[] reservedCode=new String[] {
		"%21", "%23", "%24", "%26", "%27", "%28", "%29", "%2A", "%2B", "%2C", "%2F", "%3A", "%3B", "%3D", "%3F", "%40", "%5B", "%5D"
	};
	
	public static final String[] commonSymbol=new String[] {
		"\n", " ", "\"", "%", "<", ">", "\\", "|",
	};
	
	public static final String[] commonCode=new String[] {
		"%0A", "%20", "%22", "%25", "%3C", "%3E", "%5C", "%7C"
	};

	private static final BitSet isSpecial=makeBitSet(reservedSymbol, commonSymbol);
	private static final BitSet isCommon=makeBitSet(commonSymbol);

//	public static final HashMap<String, String> codeToSymbolMap=makeCodeToSymbolMap();
//	public static final HashMap<String, String> symbolToCodeMap=makeSymbolToCodeMap();
	public static final String[] symbolToCodeArray=makeSymbolToCodeArray();
	
	/** Don't print caught exceptions */
	public static boolean suppressErrors=false;
	
}
