package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import prok.NcrnaEndpointTable;
import prok.NcrnaEndpointTable.LeaveOneOut;
import prok.TrnaBoundaryFeatures;
import prok.TrnaBoundaryFeatures.BoundaryType;
import prok.TrnaNinemerTableBuilder;
import prok.TrnaNinemerTableBuilder.LoadedTable;

/** Literal endpoint counts, repeated-key leave-one-out and real caller-table loading.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEndpointTableTest{

	public static void main(String[] args) throws Exception{
		assert(args.length==0) : "Endpoint fixture is finite and takes no external inputs";
		final byte[] repeated=bases("ACGACGACG"), other=bases("TTT"), ambiguous=bases("NAA");
		final NcrnaEndpointTable table=new NcrnaEndpointTable(BoundaryType.START, 3, 0);
		final NcrnaEndpointTable worker=new NcrnaEndpointTable(BoundaryType.START, 3, 0);
		table.add(repeated, 0);
		worker.add(repeated, 3);
		worker.add(other, 0);
		worker.add(ambiguous, 0);
		worker.add(repeated, 8);
		table.merge(worker);
		assert(table.counted()==3 && table.missing()==1 && table.ambiguous()==1)
			: "Five input endpoints comprise three valid words, one context edge and one ambiguous word";
		assert(table.rawCount(repeated, 0)==2 && table.rawCount(repeated, 3)==2 && table.rawCount(other, 0)==1)
			: "Merging adds word frequencies, not candidate positions";
		assert(table.rawCount(repeated, -1)==0 && table.rawCount(repeated, 8)==0 && table.rawCount(ambiguous, 0)==0)
			: "Invalid or ambiguous candidate words have no observed count";
		final LeaveOneOut loo=table.leaveOneOut();
		loo.set(repeated, 0);
		assert(loo.rawCount(repeated, 0)==1 && loo.rawCount(repeated, 3)==1 && loo.rawCount(repeated, 6)==1)
			: "A true ACG contribution is subtracted at every matching lookup, including positions beyond true-boundary +/-2";
		assert(loo.rawCount(other, 0)==1 && table.rawCount(repeated, 0)==2)
			: "Leave-one-out must not change another word or mutate the shared parent table";
		final float[] features=TrnaBoundaryFeatures.enrichmentProfile(repeated, 3, BoundaryType.START, 3, 0, loo);
		assert(features[0]==0.125f && features[1]==1 && features[2]==0.5f)
			: "Repeated-word candidate has counts [0,0,1,0,0], hence F=2, log2(F)/8=1/8, center=1 and rival=1/2";
		loo.set(ambiguous, 0);
		assert(loo.rawCount(repeated, 0)==2) : "An ambiguous true word contributed nothing and subtracts nothing";
		loo.set(repeated, 8);
		assert(loo.rawCount(repeated, 0)==2) : "An out-of-context true word contributed nothing and subtracts nothing";
		fails(()->loo.set(bases("CCC"), 0));

		final NcrnaEndpointTable stop=new NcrnaEndpointTable(BoundaryType.STOP, 2, 1);
		stop.add(repeated, 1);
		assert(stop.rawCount(repeated, 0)==1 && stop.rawCount(repeated, 1)==0)
			: "STOP includes the final inside base: inside2/outside1 at boundary1 counts positions0..2 (ACG)";
		final NcrnaEndpointTable start=new NcrnaEndpointTable(BoundaryType.START, 1, 2);
		start.add(repeated, 2);
		assert(start.rawCount(repeated, 0)==1) : "START inside1/outside2 at boundary2 counts positions0..2 (ACG)";
		final NcrnaEndpointTable longest=new NcrnaEndpointTable(BoundaryType.START, 15, 0);
		final byte[] allT=bases("TTTTTTTTTTTTTTT");
		longest.add(allT, 0);
		assert(longest.rawCount(allT, 0)==1) : "The largest allowed word uses thirty positive key bits without colliding with a map sentinel";

		final Path root=Files.createTempDirectory("endpoint-table-");
		try(AutoCloseable cleanup=()->{
			Files.deleteIfExists(root.resolve("start.tsv.gz"));
			Files.delete(root);
		}){
			final String path=root.resolve("start.tsv.gz").toString();
			table.write(path);
			final LoadedTable loaded=TrnaNinemerTableBuilder.loadTable(path);
			assert(loaded.k==3 && loaded.type==BoundaryType.START && loaded.insideCount==3 && loaded.outsideCount==0)
				: "The actual caller loader must recover the writer's endpoint geometry";
			assert(loaded.table.rawCount(repeated, 3)==2 && loaded.table.rawCount(other, 0)==1)
				: "The actual caller loader must recover literal raw counts from compressed output";
		}
		fails(()->new NcrnaEndpointTable(BoundaryType.START, 0, 0));
		fails(()->new NcrnaEndpointTable(BoundaryType.START, -1, 4));
		fails(()->new NcrnaEndpointTable(BoundaryType.START, Integer.MAX_VALUE, 1));
		fails(()->table.merge(stop));
		fails(()->table.merge(table));
		System.out.println("Endpoint tables passed: literal counts, worker merge, start/stop geometry, repeated-key leave-one-out, sequence limitations and native compressed loading");
	}

	private static byte[] bases(final String sequence){return sequence.getBytes(StandardCharsets.US_ASCII);}

	private static void fails(final Runnable action){
		try{
			action.run();
		}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid geometry, merge or leave-one-out contribution must fail loudly");
	}
}
