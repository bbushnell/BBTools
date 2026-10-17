package ddl;

import java.util.ArrayList;
import java.util.Arrays;
import fileIO.FileFormat;
import shared.Shared;
import stream.Read;
import stream.StreamerFactory;

/** Exercises both direct SSU consumer loaders without starting a server or loading a database.
 * @author Brian, Raiden */
public final class SSUFamilyResourceTest {
	public static void main(String[] args){
		if(args.length!=1 || !args[0].startsWith("baseline=")){
			throw new IllegalArgumentException("Supply baseline=/absolute/pre-migration/resources");
		}
		Shared.setThreads(1);Read.TO_UPPER_CASE=true;
		final String path=args[0].substring(9)+"/16S_consensus_sequence.fa";
		final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		final ArrayList<Read> reads=StreamerFactory.getReads(-1, false, ff, null, null, null);
		assert(reads!=null && reads.size()==5 && reads.get(0).id.startsWith("universal"))
			: "The baseline must be the five-model 16S library whose first reference drives SSU classification";
		final byte[] first=reads.get(0).bases;
		assert(first!=null && first.length>0 && !Arrays.equals(first, reads.get(4).bases))
			: "The fixture must distinguish the first reference from another model or an empty load";
		assert(Arrays.equals(first, SSUCompare.load16SConsensus()))
			: "SSUCompare sequence and literal classification must receive exactly the old first reference";
		assert(Arrays.equals(first, SSUServer.load16SConsensus()))
			: "SSUServer classification must receive exactly the old first reference";
		System.out.println("PASS SSU direct consumers: both loaders retain the exact first universal 16S reference");
	}
}
