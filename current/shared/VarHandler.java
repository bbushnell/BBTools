package shared;

import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;

import structures.ByteBuilder;

/**
 * VarHandle-based optimizations for Java 9+.
 * Isolated to avoid compilation issues on Java 7-8.
 * @author Brian Bushnell & Isla
 * @date January 2025
 */
final class VarHandler{
	
	private static final VarHandle STRING_VALUE_HANDLE;
	private static final VarHandle STRING_CODER_HANDLE;
	static final boolean AVAILABLE;
	
	static{
		VarHandle tempValueHandle=null;
		VarHandle tempCoderHandle=null;
		boolean success=false;
		
		try{
			MethodHandles.Lookup lookup=MethodHandles.privateLookupIn(String.class, MethodHandles.lookup());
			tempValueHandle=lookup.findVarHandle(String.class, "value", byte[].class);
			tempCoderHandle=lookup.findVarHandle(String.class, "coder", byte.class);
			success=true;
		}catch(Exception | Error e){
			// VarHandles not available - will fall back to other methods
		}
		
		STRING_VALUE_HANDLE=tempValueHandle;
		STRING_CODER_HANDLE=tempCoderHandle;
		AVAILABLE=success;
	}
	
	/**
	 * Append a String to ByteBuilder using VarHandles for direct access.
	 * Assumes ASCII/Latin-1 encoding.
	 * @param bb ByteBuilder to append to
	 * @param x String to append
	 */
	static ByteBuilder appendString(ByteBuilder bb, String x){
		if(x==null){return bb;}
		
		byte[] value=(byte[])STRING_VALUE_HANDLE.get(x);
		byte coder=(byte)STRING_CODER_HANDLE.get(x);
		
		if(coder==0){ // LATIN1 - direct copy
			bb.expand(value.length);
			System.arraycopy(value, 0, bb.array, bb.length, value.length);
			bb.length+=value.length;
		}else{ // UTF-16 fallback
			bb.expand(x.length());
			for(int i=0; i<x.length(); i++){
				bb.array[bb.length++]=(byte)x.charAt(i);
			}
		}
		return bb;
	}
	
	private VarHandler(){} // Prevent instantiation
}