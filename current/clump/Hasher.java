package clump;

import java.util.Random;

import shared.Tools;
import stream.Read;

public class Hasher {

	private static synchronized long[][] makeCodes2(int modes){
		long[][] r=makeCodes(128, modes);
		
		for(int i=0; i<26; i++){
			char c=(char)('A'+i);
			r[Tools.toLowerCase(c)]=r[c];
		}
		return r;
	}
	
	private static synchronized long[][] makeCodes(int symbols, int modes){
		Random randy=new Random(1);
		long[][] r=new long[symbols][modes];
		for(int i=0; i<symbols; i++){
			for(int j=0; j<modes; j++){
				r[i][j]=randy.nextLong();
			}
		}
		return r;
	}
	
	public static long hash(byte[] bases){
		long code=bases.length;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			int mode=(int)(code&31);
			assert(hashcodes[b]!=null) : "Invalid sequence character: '"+(char)b+"'";
			code=code^hashcodes[b][mode];
			code=Long.rotateLeft(code, 1);
		}
		return code;
	}
	
	public static final long hash(Read r){
		return hash(r.bases);
	}
	
	public static final long hashPair(Read r){
		long a=hash(r);
		if(r.mate==null){return a;}
		long b=hash(r.mate);
		return a^Long.rotateLeft(b, 1);
	}
	
	public static final boolean equalsPaired(Read a, Read b){
		return equals(a, b) && equals(a.mate, b.mate);
	}
	
	public static final boolean equals(Read a, Read b){
		if(a==b){return true;}
		if(a==null || b==null){
			assert(a!=null || b!=null);
			return false;
		}
		if(a.length()!=b.length()){return false;}
		if(a.length()==0){return true;}
		byte[] ab=a.bases, bb=b.bases;
		for(int i=0; i<ab.length; i++){
			if(ab[i]!=bb[i]){return false;}
		}
		return true;
	}
	
	private static final long[][] hashcodes=makeCodes2(32);
	
}
