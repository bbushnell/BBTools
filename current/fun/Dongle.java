package fun;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

import shared.Tools;

public class Dongle {
	
	public static void main(String[] args) {
		long millis=System.currentTimeMillis();
		if(args.length>0) {
			millis=parse(args[0]);
		}
		long a=encode(millis);
		System.err.println(a);
		System.err.println(millis);
		System.err.println(decode(a));
	}
	
	public static final boolean check(Object...args) {
		long a=min, b=max;
		if(args!=null && args.length>1) {a=(Long)args[1];}
		if(args!=null && args.length>2) {b=(Long)args[2];}
		a=decode(a);
		b=decode(b);
		long x=System.currentTimeMillis();
		return (b>a && b-a<limit && x>=a && x<=b);
	}
	
	private static long encode(long x) {
		x^=number;
		long a=x&mask;
		long b=x&(mask<<1);
		x=a|(Long.rotateLeft(b, rot));
		a=x&mask2;
		b=x&(mask2<<2);
		x=a|Long.rotateRight(b, rot2);
		a=x&mask3;
		b=x&(mask3>>>4);
		x=a|Long.rotateLeft(b, rot3);
		return x;
	}
	
	private static long decode(long x) {
		long a=x&mask3;
		long b=x&(mask3>>>4);
		x=a|(Long.rotateRight(b, rot3));
		a=x&mask2;
		b=x&(mask2<<2);
		x=(a|(Long.rotateLeft(b, rot2)));
		a=x&mask;
		b=x&(mask<<1);
		x=(a|(Long.rotateRight(b, rot)))^number;
		return x;
	}
	
	private static long parse(String s) {
		try {
			if(Tools.isNumeric(s)) {return Long.parseLong(s);}
			Date d=sdf.parse(s);
			System.err.println(d);
			return d.getTime();
		} catch (ParseException e) {
			return -1L;
		}
	}
	
	private static final String pattern="MM-dd-yyyy";
	private static final SimpleDateFormat sdf=new SimpleDateFormat(pattern);
	private static final long min=-6788251374689658131L;
	private static final long max=2715152938288332905L;
	private static final long number=4964420948893066024L;
	private static final long mask=0x5555555555555555L;
	private static final long mask2=0x3333333333333333L;
	private static final long mask3=0xF0F0F0F0F0F0F0F0L;
	private static final long limit=346896001029L;
	private static final int rot=26;
	private static final int rot2=44;
	private static final int rot3=16;
	
}
