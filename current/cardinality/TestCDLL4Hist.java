package cardinality;
public class TestCDLL4Hist {
	public static void main(String[] args){
		CompressedDynamicLogLog4.DISABLE_EEMASK=false;
		if(args.length>0 && args[0].equals("noeemask")){CompressedDynamicLogLog4.DISABLE_EEMASK=true;}
		CompressedDynamicLogLog4.DEBUG_HIST=false;
		final CompressedDynamicLogLog4 c=new CompressedDynamicLogLog4(512, 31, 12345L, 0);
		final java.util.Random r=new java.util.Random(42);
		for(int i=0; i<5000000; i++){
			c.add(r.nextLong());
		}
		double[] est=c.rawEstimates();
		System.err.println("minZeros="+c.getMinZeros()+" card="+c.cardinality());
	}
}
