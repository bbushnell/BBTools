package test.prok;


import java.util.Arrays;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import prok.RrnaResourceIO;
import structures.ByteBuilder;

/** Saved two-radius candidate blocks retain zero-positive ends and edge clipping.
 * Preserved mitochondrial endpoint experiment: this class supports construction,
 * inference or validation of the experimental endpoint method. The method was not
 * promoted, and the mitochondrial-SSU direction was set aside on October 24 because
 * these genes did not fit the prokaryotic calling approach. Kept for a future
 * mitochondrial calling investigation; ordinary CallGenes never enables it.
 * @author Brian Bushnell, Raiden */
public final class MitoEndpointTrainingVectorsTest {
	public static void main(String[] args){
		check(args.length==1,"Provide a fresh output prefix within the run directory");
		final String small=args[0]+".r128.tsv",large=args[0]+".r256.tsv";
		check(!new java.io.File(small).exists() && !new java.io.File(large).exists(),"The vector fixture must not overwrite existing evidence");
		final byte[] seq=new byte[1024];Arrays.fill(seq,(byte)'A');
		final MitoEndpointTables tables=new MitoEndpointTables(1,9,9,1000);
		tables.add(seq,450,900,0);final long[][] counts=new long[2][4],curve=new long[7][4];
		final ByteStreamWriter a=RrnaResourceIO.writer(small),b=RrnaResourceIO.writer(large);
		try{MitoEndpointTrainingVectors.writeCandidates(seq,300,800,450,900,0,.5f,.75f,false,tables,a,b,new ByteBuilder(),counts,curve);}
		finally{RrnaResourceIO.finish(a,b);}
		check(Arrays.equals(counts[0],new long[]{257,257,0,1}),"Radius128 must retain the start's zero-positive block and the stop's positive block");
		check(Arrays.equals(counts[1],new long[]{513,480,1,1}),"Radius256 must reach both truths and clip stop candidates at the actual window edge");
		check(Arrays.equals(read(small),counts[0]) && Arrays.equals(read(large),counts[1]),"Saved rows must conserve both radius counters despite reusing the output builder");
		check(Arrays.equals(curve[4],counts[0]) && Arrays.equals(curve[6],counts[1]),"Curve counters must equal actual generated128/256datasets");
		check(Arrays.equals(curve[0],new long[]{65,65,0,0}) && Arrays.equals(curve[5],new long[]{385,385,1,1}),"Intermediate32/192radius reachability must follow actual candidate positions, not interpolate the endpoint curve");
		System.out.println("PASS MitoEndpointTrainingVectorsTest zero_positive_blocks_two_radii_edge_clipping_saved_rows");
	}
	private static long[] read(String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final long[] n=new long[4];
		try{for(byte[] line;(line=in.nextLine())!=null;){p.set(line);check(p.terms()==11,"Each saved candidate needs ten features plus label");
			final float flag=p.parseFloat(3);check(flag==0 || flag==1,"V3 endpoint flag must preserve its column");final int end=(int)flag,label=p.parseInt(10);
			check((label==0 || label==1) && p.parseFloat(5)==.5f && p.parseFloat(6)==.75f,"Source GC, genomic identity and binary target must survive serialization");n[end]++;n[end+2]+=label;}}
		finally{check(!in.close(),"Saved vector fixture read failed");}return n;
	}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
