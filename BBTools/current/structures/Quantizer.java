package structures;

import java.util.ArrayList;

import shared.KillSwitch;
import shared.Parse;
import shared.Tools;
import stream.Read;

public class Quantizer {
	
	public static boolean parse(String arg, String a, String b){
		if(a.equals("quantize")){
			if(b!=null && b.equalsIgnoreCase("sticky")){
				STICKY=true;
				return true;
			}
		}else if(a.equals("quantizesticky")){
			if(b!=null && (b.charAt(0)=='/' || Tools.isDigit(b.charAt(0)))) {
				STICKY=true;
			}else {
				STICKY=Parse.parseBoolean(b);
				return true;
			}
		}
		
		if(b==null || b.length()<1 || Character.isLetter(b.charAt(0))){
			return Parse.parseBoolean(b);
		}
		return setArray(b);
	}
	
	private static boolean setArray(String s){
		final byte[] array;
		if(s.indexOf(',')<0){
			if(s.charAt(0)=='/') {s=s.substring(1);}
			int quant=Integer.parseInt(s);
			assert(quant>0 && quant<128);
			if(quant==1){return false;}
			ByteBuilder bb=new ByteBuilder();
			for(int i=0, max=Read.MAX_CALLED_QUALITY(); i<=max; i+=quant){
				bb.append((byte)i);
			}
			array=bb.toBytes();
		}else{
			array=Parse.parseByteArray(s, ",");
		}
		setArray(array);
		return true;
	}
	
	private static void setArray(byte[] a){
		quantizeArray=a;
		qualityRemapArray=makeQualityRemapArray(quantizeArray);
	}
	
	public static void quantize(ArrayList<Read> list) {
		if(list==null) {return;}
		for(Read r : list) {
			if(r!=null) {
				quantize(r);
				quantize(r.mate);
			}
		}
	}
	
	public static void quantize(Read r1, Read r2){
		quantize(r1);
		quantize(r2);
	}
	
	public static void quantize(Read r){
		if(r!=null) {quantize(r.quality);}
	}
	
	public static void quantize(byte[] quals){
		if(quals==null){return;}
		byte prev=0;
		for(int i=0; i<quals.length; i++){
			final byte qOld=quals[i];
			final byte q0=qualityRemapArray[qOld];
			byte q=q0;
			if(STICKY && q!=prev && prev>0 && q>0 && Tools.absdif(qOld, prev)<=Tools.absdif(qOld, q)){q=prev;}
//			assert(q==q0) : STICKY+", "+qOld+" -> "+q0+" -> "+q+", prev="+prev;
			quals[i]=q;
			prev=q;
		}
	}
	
	private static final byte[] makeQualityRemapArray(byte[] quantizeArray) {
		byte[] array=KillSwitch.allocByte1D(128);
		for(int i=0; i<array.length; i++){
			byte q=0;
			for(byte x : quantizeArray){
				if((i>0 && q==0 && x>0) || Tools.absdif(x, i)<=Tools.absdif(q, i)){q=x;}
			}
			array[i]=q;
		}
		return array;
	}
	
//	private static byte[] quantizeArray={0, 8, 13, 22, 27, 32, 37}; //Old
	private static byte[] quantizeArray={0, 14, 21, 27, 32, 36};
	private static byte[] qualityRemapArray=makeQualityRemapArray(quantizeArray);
	private static boolean STICKY=true;
	
}
