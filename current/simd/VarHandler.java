package simd;

import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;

import structures.ByteBuilder;

/**
 * Optional Java9+ access to String storage for ByteBuilder appends.
 * Vector isolates this backend from older JVMs. Access also requires a String
 * implementation with byte[] value and byte coder fields, and permission to
 * obtain a private lookup into java.lang. Failure leaves AVAILABLE false.
 * String storage is borrowed read-only; the caller owns the destination builder.
 * @author Brian Bushnell, Isla
 * @date January 2025
 */
final class VarHandler{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Prevents construction of this static backend. */
	private VarHandler(){}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Appends the low byte of each UTF-16 code unit, matching the scalar writer.
	 * Latin-1 storage is copied directly; other strings use charAt and byte casts.
	 * This does not encode arbitrary text as UTF-8. AVAILABLE must be true before
	 * calling with non-null text. The destination must not be shared concurrently.
	 * @param bb ByteBuilder to append to
	 * @param x String to append; null leaves the destination unchanged
	 * @return The same destination builder
	 */
	static ByteBuilder appendString(ByteBuilder bb, String x){
		if(x==null){return bb;}

		//String owns this array. Read it only; mutating it would violate String immutability.
		//Older JVMs must avoid linking this backend; catching lookup failure cannot supply a missing VarHandle API.
		byte[] value=(byte[])STRING_VALUE_HANDLE.get(x);
		byte coder=(byte)STRING_CODER_HANDLE.get(x);

		if(coder==0){//Latin-1: one byte per code unit
			bb.expand(value.length);
			System.arraycopy(value, 0, bb.array, bb.length, value.length);
			bb.length+=value.length;
		}else{//UTF-16: retain the scalar writer's low-byte conversion
			bb.expand(x.length());
			for(int i=0; i<x.length(); i++){
				bb.array[bb.length++]=(byte)x.charAt(i);
			}
		}
		return bb;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	/** Read-only access to String's backing bytes, or null if lookup failed. */
	private static final VarHandle STRING_VALUE_HANDLE;
	/** Read-only access to String's Latin-1/UTF-16 selector, or null on failure. */
	private static final VarHandle STRING_CODER_HANDLE;
	/** True only when both private field lookups succeeded. */
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
		}catch(Exception|Error e){
			// VarHandles not available - will fall back to other methods
		}

		STRING_VALUE_HANDLE=tempValueHandle;
		STRING_CODER_HANDLE=tempCoderHandle;
		AVAILABLE=success;
	}

}
