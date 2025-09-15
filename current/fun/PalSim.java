package fun;

import java.util.Random;

public class PalSim {
	
	public static void main(String[] args) {
		float a0=Float.parseFloat(args[0]);
		float b0=Float.parseFloat(args[1]);
		float c0=Float.parseFloat(args[2]);
		if(a0>.3 || b0>.3 || c0>.3) {a0*=0.3; b0*=0.3; c0*=0.3;}
		assert(a0>=0 && a0<=.3);
		assert(b0>=0 && b0<=.3);
		assert(c0>=0 && c0<=.3);
		final long sims=(args.length>3 ? Long.parseLong(args[3]) : 1000000);
		final Random randy=new Random();
		final float product0=(1+a0)*(1+b0)*(1+c0);
		final float maxProduct=1.3f*1.3f*1.3f;
		
		long betterCount=0;
		for(long i=0; i<sims; i++) {
			final float a=randy.nextFloat()*0.3f+1f;
			final float b=randy.nextFloat()*0.3f+1f;
			final float c=randy.nextFloat()*0.3f+1f;
			final float product=a*b*c;
			betterCount+=(product0>=product ? 1 : 0);
		}
		
		final float ratio=betterCount/(float)sims;
		final float rawPlacement=(product0-1)/(maxProduct-1);
		final double third=1.0/3.0;
		final double rootPlacement=(Math.pow(product0, third)-1)/(Math.pow(maxProduct, third)-1);
		final double linearPlacement=(a0+b0+c0)/(0.9);

		System.out.println("sim:\t"+String.format("%.3f%%", ratio*100));
		System.out.println("raw:\t"+String.format("%.3f%%", rawPlacement*100));
		System.out.println("root:\t"+String.format("%.3f%%", rootPlacement*100));
		System.out.println("linear:\t"+String.format("%.3f%%", linearPlacement*100));
	}
	
}
