package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Objects;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/** Integration oracle: the legacy union and exact ordered partition packs must
 * produce identical assignments, scores, coverage and candidate accounting.
 * Uses real package readers, shortlist scoring and positional alignments.
 * @author Keqing
 */
public final class ProteinFamilyCompetitionTest {

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, ProteinFamilyCompetitionTest.class, false);
		try{run(HmmComparisonData.options(pp.args));}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final HashMap<String,String> options) throws Exception{
		for(String key : options.keySet()){
			check(Arrays.asList("models", "packs", "in", "out", "count").contains(key), "Unknown competition-test argument: "+key);
		}
		final Path out=Paths.get(HmmComparisonData.required(options, "out"));
		check(!Files.exists(out), "Competition test output must be fresh");
		final int expected=Integer.parseInt(HmmComparisonData.required(options, "count"));
		check(expected>0, "Expected query count must be positive");
		Shared.setThreads(1); Shared.AMINO_IN=true;
		final ProteinFamilyAssigner union=ProteinFamilyAssigner.load(Paths.get(HmmComparisonData.required(options, "models")));
		final String[] paths=HmmComparisonData.required(options, "packs").split(",", -1);
		check(paths.length>=2, "Partition parity requires at least two independently loaded packs");
		final ProteinFamilyAssigner[] parts=new ProteinFamilyAssigner[paths.length];
		int offset=0;
		for(int p=0; p<paths.length; p++){
			parts[p]=ProteinFamilyAssigner.load(Paths.get(paths[p]));
			for(int r=0; r<parts[p].familyCount(); r++){
				check(offset<union.familyCount() && parts[p].familyId(r)==union.familyId(offset)
					&& parts[p].representative(r).equals(union.representative(offset))
					&& parts[p].memberCount(r)==union.memberCount(offset),
					"Partitions do not preserve the union roster at rank "+offset);
				offset++;
			}
		}
		check(offset==union.familyCount(), "Partition roster omits union families");
		final ProteinFamilyCompetition combined=new ProteinFamilyCompetition(parts);
		check(combined.familyCount()==union.familyCount() && combined.packCount()==parts.length, "Combined dimensions changed");
		reject(()->new ProteinFamilyCompetition(union, union), "Duplicate representative");
		final ProteinFamilyAssigner.Worker legacy=union.newWorker();
		final ProteinFamilyCompetition.Worker single=new ProteinFamilyCompetition(union).newWorker();
		final ProteinFamilyCompetition.Worker multiple=combined.newWorker();
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(
			HmmComparisonData.required(options, "in"), FileFormat.FASTA, null, true, false), null, true, -1);
		int queries=0, assigned=0;
		input.start();
		try{
			for(ListNum<Read> batch=input.nextList(); batch!=null && batch.size()>0; batch=input.nextList()){
				for(Read read : batch){
					check(read.mate==null, "Competition fixture expects unpaired protein sequences");
					final String id=HbmCompetitiveAssign.queryId(read.id);
					final ProteinFamilyAssigner.Hit a=legacy.assign(id, read.bases);
					final ProteinFamilyAssigner.Hit b=single.assign(id, read.bases);
					final ProteinFamilyAssigner.Hit c=multiple.assign(id, read.bases);
					equal(a, b, id+" single"); equal(a, c, id+" partitioned");
					if(a.assigned()){
						check(b.packIndex==0 && b.rank==a.rank, "Single-pack rank mapping changed for "+id);
						check(c.packIndex>=0 && c.packIndex<parts.length, "Missing winner pack for "+id);
						int rank=c.rank;
						for(int p=0; p<c.packIndex; p++){rank+=parts[p].familyCount();}
						check(rank==a.rank && parts[c.packIndex].familyId(c.rank)==c.familyId, "Wrong local winner rank for "+id);
						assigned++;
					}else{check(b.packIndex==-1 && c.packIndex==-1, "No-hit result retained a winning pack for "+id);}
					queries++;
				}
			}
		}finally{input.close();}
		check(!input.errorState() && queries==expected, "Query reader failed or count differs: "+queries+" != "+expected);
		final ByteStreamWriter writer=HmmComparisonData.writer(out.toString());
		writer.println("FAMILY_COMPETITION_PARITY_PASS queries="+queries+" assigned="+assigned
			+" packs="+parts.length+" legacy_single_exact=true partition_exact=true duplicate_rejected=true");
		HmmComparisonData.close(writer);
		System.err.println("FAMILY_COMPETITION_PARITY_PASS queries="+queries+" assigned="+assigned+" packs="+parts.length);
	}

	private static void equal(final ProteinFamilyAssigner.Hit a, final ProteinFamilyAssigner.Hit b, final String id){
		check(a.assigned()==b.assigned() && a.familyId==b.familyId && Objects.equals(a.representative, b.representative)
			&& a.queryLength==b.queryLength && a.score64==b.score64 && a.referenceStart==b.referenceStart
			&& a.referenceEnd==b.referenceEnd && Float.compare(a.identity, b.identity)==0
			&& Float.compare(a.queryCoverage, b.queryCoverage)==0 && Float.compare(a.coreCoverage, b.coreCoverage)==0
			&& Float.compare(a.mutualCoverage, b.mutualCoverage)==0 && a.scored==b.scored
			&& a.lengthRejected==b.lengthRejected && b.scored+b.lengthRejected==50,
			"Combined competition differs from legacy assignment for "+id);
	}

	private static void reject(final Runnable action, final String reason){
		try{action.run();}catch(IllegalArgumentException e){check(e.getMessage().contains(reason), "Unexpected rejection: "+e); return;}
		throw new AssertionError("Expected rejection: "+reason);
	}
	private static void check(final boolean ok, final String reason){if(!ok){throw new AssertionError(reason);}}
}
