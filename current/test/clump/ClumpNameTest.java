package test.clump;

import java.nio.charset.StandardCharsets;

import clump.Clumpify;
import stream.Read;

/** Checks numeric short names and the nullable-mate calling convention.
 * @author Shinobu
 */
public class ClumpNameTest{

	/** Runs baseline|fixed null handling plus unchanged literal name controls. */
	public static void main(String[] args){
		assert(args.length==1) : "Expected baseline or fixed";
		final boolean fixed=args[0].equals("fixed");
		assert(fixed || args[0].equals("baseline")) : "Unknown expectation: "+args[0];
		boolean rejected=false;
		try{Clumpify.shortName(null);}catch(NullPointerException e){rejected=true;}
		assert(rejected!=fixed) : "Unpaired callers pass r.mate==null to shortName";
		System.out.println("case\tvalue");
		System.out.println("null-mate\t"+(rejected ? "NullPointerException" : "no-op"));
		final long[] ids={0, 1, 9, 10, 999999999L, 1000000000L, Long.MAX_VALUE};
		final String[] names={"0000000000", "0000000001", "0000000009", "0000000010",
				"0999999999", "1000000000", "9223372036854775807"};
		for(int i=0; i<ids.length; i++){
			for(int pair=0; pair<2; pair++){
				final Read read=new Read("ACGT".getBytes(StandardCharsets.US_ASCII), null, "old", ids[i]);
				read.setPairnum(pair);
				Clumpify.shortName(read);
				assert(read.id.equals(names[i]+" "+(pair+1)+":")) : "Wrong short name: "+read.id;
				System.out.println(ids[i]+"-mate"+pair+"\t"+read.id);
			}
		}
		Clumpify.shrinkName(null);
		final Read read=new Read("ACGT".getBytes(StandardCharsets.US_ASCII), null,
				"HISEQ:000000000-MISEQ:1:2:3", 0);
		Clumpify.shrinkName(read);
		assert(read.id.equals("H:M:1:2:3")) : "Sequencer-name shrink control changed: "+read.id;
		System.out.println("shrink-control\t"+read.id);
	}
}
