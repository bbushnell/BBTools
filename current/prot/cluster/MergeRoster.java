package prot.cluster;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import fileIO.ByteFile;
import map.LongHashSet;
import parse.LineParser1;
import prot.ProteinSequence;

/** Original-tail weights bound to the ordered Stage-B consensus inputs.
 * Reads the roster emitted by DdlClustStageB.collect, never a guessed FASTA size.
 * @author Keqing
 */
public final class MergeRoster{

	public static int[] load(String path,ProteinSequence[] sequences){
		if(path==null || sequences==null){throw new IllegalArgumentException("Merge roster needs a path and ordered consensus inputs");}
		final int[] weights=new int[sequences.length];
		final LongHashSet ids=new LongHashSet();
		final LineParser1 parser=new LineParser1('\t');
		final ByteFile file=ByteFile.makeByteFile(path,false);
		int ordinal=0;
		try{
			if(!Arrays.equals(file.nextLine(),HEADER)){
				throw new IllegalArgumentException("Unexpected Stage-B merge roster header: "+path);
			}
			for(byte[] line=file.nextLine(); line!=null; line=file.nextLine(), ordinal++){
				parser.set(line);
				if(parser.terms()!=5 || ordinal>=weights.length || nonnegative(parser,0)!=ordinal){
					throw new IllegalArgumentException("Merge roster must contain exactly one ordered row per consensus; row="+ordinal);
				}
				final int shard=nonnegative(parser,2), family=nonnegative(parser,3), count=nonnegative(parser,4);
				final long key=((long)shard<<32)|family;
				final String expected="s"+shard+"_f"+family;
				if(count<1 || sequences[ordinal]==null || !expected.equals(sequences[ordinal].id) ||
					!parser.termEquals(expected,1) || !ids.add(key)){
					throw new IllegalArgumentException("Merge roster needs unique matching family IDs and positive original-tail counts; row="+ordinal);
				}
				weights[ordinal]=count;
			}
		}finally{
			if(file.close()){throw new IllegalStateException("Failed reading merge roster: "+path);}
		}
		if(ordinal!=weights.length){throw new IllegalArgumentException("Merge roster omitted consensus inputs: "+ordinal+" != "+weights.length);}
		assert(ordinal==sequences.length) : "Each merge query must retain its original-tail weight at the same input ordinal";
		return weights;
	}

	/** Checked decimal parsing prevents overflow from turning a large family into a small one. */
	private static int nonnegative(LineParser1 parser,int term){
		parser.setBounds(term);
		final byte[] line=parser.line();
		final int start=parser.a(), end=parser.b();
		if(start>=end){throw new IllegalArgumentException("Empty numeric merge roster field "+term);}
		int value=0;
		for(int i=start; i<end; i++){
			final int digit=line[i]-'0';
			if(digit<0 || digit>9 || value>(Integer.MAX_VALUE-digit)/10){
				throw new IllegalArgumentException("Invalid or overflowing numeric merge roster field "+term);
			}
			value=value*10+digit;
		}
		assert(value>=0) : "Merge roster sizes and canonical source IDs are nonnegative signed integers";
		return value;
	}

	private static final byte[] HEADER="stagea_ordinal\tstagea_id\tshard\tstagea_family_id\tsource_members".getBytes(StandardCharsets.US_ASCII);
}
