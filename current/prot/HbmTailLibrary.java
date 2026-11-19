package prot;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;
import template.ThreadWaiter;

/** Prepares active stop-tail singletons for the existing covering-set and F4 index builders.
 * Worker threads verify and copy independent families. Retired groups remain in summary denominators;
 * original unclustered/orphan queries remain in the separate full build-tail FASTA.
 * Source singleton graphs are verified and never rewritten.
 * @author Keqing
 */
public final class HbmTailLibrary {

	public static void main(String[] args){
		try{
			final PreParser pp=new PreParser(args, HbmTailLibrary.class, false);
			try{runForTest(pp.args);}finally{Shared.closeStream(pp.outstream);}
		}catch(Throwable failure){failure.printStackTrace(); System.exit(1);}
	}

	/** Throwing entry point also permits malformed-input fixtures to verify absence of PASS. */
	public static void runForTest(String[] args) throws Exception{
		Shared.AMINO_IN=true;
		final HashMap<String,String> options=HmmComparisonData.options(args);
		for(String key : options.keySet()){
			require(Arrays.asList("families", "familiessha80", "background", "backgroundsha80", "out", "t").contains(key), "Unknown tail-library option: "+key);
		}
		final int threads=Integer.parseInt(options.getOrDefault("t", "1"));
		require(threads>0, "Tail-library threads must be positive: "+threads);
		Shared.setThreads(threads);
		final String input=HmmComparisonData.required(options, "families"), inputPin=HmmComparisonData.required(options, "familiessha80");
		final String bg=HmmComparisonData.required(options, "background"), bgPin=HmmComparisonData.required(options, "backgroundsha80");
		HbmProfilePilot.requireHash(input, inputPin); HbmProfilePilot.requireHash(bg, bgPin);
		NewFamilyHbmBuilder.readBackground(bg);
		final List<String[]> rows=HmmComparisonData.rows(input);
		require(rows.size()>1 && String.join("\t", rows.get(0)).equals(REGISTRY_HEADER), "Expected the tail-builder family registry");
		final Path out=Paths.get(HmmComparisonData.required(options, "out")).toAbsolutePath();
		require(!Files.exists(out), "Tail library output must be fresh: "+out);
		final HashSet<String> seen=new HashSet<String>();
		long proposed=0, retained=0; int active=0;
		final ArrayList<Job> jobs=new ArrayList<Job>();
		for(int i=1; i<rows.size(); i++){
			final String[] row=rows.get(i);
			require(row.length==6 && token(row[0]) && seen.add(row[0]), "Malformed or duplicate tail-family registry row: "+i);
			final long count=Long.parseLong(row[1]), kept=Long.parseLong(row[2]), rejected=Long.parseLong(row[3]);
			// An outer round may assign no copies to a formerly active family; keep its retirement row.
			require(count>=0 && kept>=0 && kept<=count && rejected==count-kept, "Tail family counts do not conserve input copies: "+row[0]);
			require(Paths.get(row[5]).isAbsolute(), "Tail family directory must be absolute: "+row[5]);
			if(row[4].equals("ACTIVE")){
				require(kept>=3, "Active tail family has fewer than three copies: "+row[0]);
				jobs.add(new Job(active, row[0], row[2], Paths.get(row[5]).toRealPath())); active++;
			}
			else{require(kept==0 && !row[4].isEmpty(), "Retired tail family cannot retain copies: "+row[0]);}
			proposed=Math.addExact(proposed, count); retained=Math.addExact(retained, kept);
		}
		require(active>0, "No active tail families; retain the seed census and diagnose before indexing");
		Files.createDirectory(out); Files.createDirectory(out.resolve("members"));
		final Result[] results=export(out, jobs, bgPin);
		final ByteStreamWriter fasta=HmmComparisonData.writer(out.resolve("consensus.faa").toString());
		final ByteStreamWriter models=HmmComparisonData.writer(out.resolve("models.tsv").toString());
		final ByteStreamWriter families=HmmComparisonData.writer(out.resolve("covering_families.tsv").toString());
		models.println(MODELS_HEADER);
		long residues=0;
		try{
			for(int rank=0; rank<results.length; rank++){
				final Result result=results[rank];
				require(result!=null && result.rank==rank, "Dense index must contain every active tail-family rank: "+rank);
				fasta.print(result.consensus);
				families.println(result.family);
				models.println(result.model);
				residues=Math.addExact(residues, result.residues);
			}
		}finally{HmmComparisonData.close(fasta); HmmComparisonData.close(models); HmmComparisonData.close(families);}
		require(results.length==active, "Dense index must contain exactly the active tail roster");
		HbmProfilePilot.requireHash(input, inputPin); HbmProfilePilot.requireHash(bg, bgPin);
		Files.copy(Paths.get(bg), out.resolve("background.tsv"));
		HbmProfilePilot.requireHash(out.resolve("background.tsv").toString(), bgPin);
		final ByteStreamWriter summary=HmmComparisonData.writer(out.resolve("summary.tsv").toString());
		summary.println("proposed_families\tactive_families\tretired_families\tproposed_copies\tretained_copies\trejected_copies\tretained_residues");
		summary.println(new ByteBuilder().append(rows.size()-1).tab().append(active).tab().append(rows.size()-1-active).tab().append(proposed)
			.tab().append(retained).tab().append(proposed-retained).tab().append(residues)); HmmComparisonData.close(summary);
		final ByteStreamWriter recipe=HmmComparisonData.writer(out.resolve("recipe.tsv").toString());
		recipe.println("families_sha80\t"+inputPin+"\nbackground_sha80\t"+bgPin+"\noperation\tverified_singletons_unchanged\nmember_policy\taccepted_original_tails");
		for(String file : new String[]{"models.tsv", "consensus.faa", "covering_families.tsv"}){recipe.println(file+"_sha80\t"+DigestSuffix.file(out.resolve(file).toString()));}
		HmmComparisonData.close(recipe);
		final ByteStreamWriter pass=HmmComparisonData.writer(out.resolve("PASS").toString());
		pass.println("HBM_TAIL_LIBRARY_PASS"); HmmComparisonData.close(pass);
	}

	private static Result[] export(final Path out, final ArrayList<Job> jobs, final String bgPin){
		final Result[] results=new Result[jobs.size()];
		final AtomicInteger next=new AtomicInteger();
		final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
		final int threads=Math.min(Shared.threads(), jobs.size());
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){workers.add(new ProcessThread(out, jobs, results, bgPin, next, failure));}
		ThreadWaiter.startAndWait(workers);
		if(failure.get()!=null){throw new RuntimeException("Tail-library worker failed", failure.get());}
		return results;
	}

	private static Result exportJob(final Path out, final Job job, final String bgPin) throws Exception{
		final Path directory=job.directory;
		final String[] pins=new String[ARTIFACTS.length];
		for(int j=0; j<pins.length; j++){pins[j]=DigestSuffix.file(directory.resolve(ARTIFACTS[j]).toString());}
		HbmProfilePilot.requireHash(directory.resolve("background.tsv").toString(), bgPin);
		final List<ProteinSequence> refs=readTails(directory.resolve("consensus.faa").toString());
		require(refs.size()==1 && refs.get(0).id.equals(job.family), "Singleton consensus does not match family: "+job.family);
		final ProteinSequence ref=refs.get(0); validateTail(ref);
		final List<String[]> roster=HmmComparisonData.rows(directory.resolve("roster.tsv").toString());
		require(roster.size()==2 && String.join("\t", roster.get(0)).equals("rank\trep_id\tmembers")
			&& Arrays.equals(roster.get(1), new String[]{"0", job.family, job.retained}), "Singleton roster differs from retained membership: "+job.family);
		final HbmBundleLoader.Loaded loaded=HbmBundleLoader.load(directory.resolve("refined.mqhb"), Collections.singletonList(job.family), id->ref.enc,
			HbmBundleLoader.loadSemanticProvenance(directory.resolve("provenance.tsv").toString()));
		require(loaded.familyCount()==1, "Expected a verified singleton graph: "+job.family);
		final List<ProteinSequence> members=readTails(directory.resolve("retained.faa").toString());
		require(members.size()==Long.parseLong(job.retained), "Accepted-member count differs from registry: "+job.family);
		final HashSet<String> memberIds=new HashSet<String>(); long familyResidues=0;
		for(ProteinSequence member : members){
			validateTail(member); require(memberIds.add(member.id), "Duplicate accepted copy: "+member.id);
			familyResidues=Math.addExact(familyResidues, member.length());
		}
		final Path memberCopy=out.resolve("members").resolve(job.rank+".faa");
		Files.copy(directory.resolve("retained.faa"), memberCopy);
		HbmProfilePilot.requireHash(memberCopy.toString(), pins[3]);
		for(int j=0; j<pins.length; j++){HbmProfilePilot.requireHash(directory.resolve(ARTIFACTS[j]).toString(), pins[j]);}
		final ByteBuilder consensus=new ByteBuilder();
		consensus.append('>').append(job.family).append(" rank=").append(job.rank).nl().append(HmmComparisonData.ascii(ref.enc)).nl();
		final ByteBuilder family=new ByteBuilder().append(job.family).tab().append(memberCopy.toString());
		final ByteBuilder model=new ByteBuilder().append(job.rank).tab().append(job.family).tab().append(directory.toString()).tab().append(job.retained).tab().append(familyResidues);
		for(String pin : pins){model.tab().append(pin);}
		return new Result(job.rank, familyResidues, consensus, family, model);
	}

	/** Strict family reader: no skipped stop-bearing records or stripped terminal markers. */
	static List<ProteinSequence> readTails(String file){
		final ArrayList<ProteinSequence> result=new ArrayList<ProteinSequence>();
		final Streamer input=StreamerFactory.makeStreamer(FileFormat.testInput(file, FileFormat.FASTA, null, true, false), 0, true, -1, false, true, 0);
		input.start();
		try{
			for(ListNum<Read> list=input.nextList(); list!=null && list.size()>0; list=input.nextList()){
				for(Read read : list){
					require(read.mate==null && read.bases!=null && read.bases.length>0, "Tail library requires nonempty unpaired protein records");
					final ProteinSequence p=new ProteinSequence(HbmCompetitiveAssign.queryId(read.id), read.bases);
					require(p.length()==read.bases.length, "Tail library cannot strip terminal stop markers: "+p.id);
					validateTail(p); result.add(p);
				}
			}
		}finally{input.close(); require(!input.errorState(), "Tail family reader failed: "+file);}
		return result;
	}

	private static void validateTail(ProteinSequence p){
		require(token(p.id) && p.length()>=100 && p.length()<=200, "Tail library needs unique-token100..200-AA sequences: "+p.id);
		for(byte b : p.enc){require(b>=0 && b<20, "Tail library excludes ambiguous or stop residues: "+p.id);}
	}
	private static boolean token(String value){
		if(value==null || value.isEmpty()){return false;}
		for(int i=0; i<value.length(); i++){if(value.charAt(i)<=32 || value.charAt(i)>=127){return false;}}return true;
	}
	private static void require(boolean ok,String message){if(!ok){throw new IllegalArgumentException(message);}}

	private static final class ProcessThread extends Thread{
		ProcessThread(final Path out, final ArrayList<Job> jobs, final Result[] results, final String bgPin,
				final AtomicInteger next, final AtomicReference<Throwable> failure){
			this.out=out; this.jobs=jobs; this.results=results; this.bgPin=bgPin; this.next=next; this.failure=failure;
			setDaemon(true);
		}
		@Override
		public void run(){
			try{
				for(int idx=next.getAndIncrement(); idx<jobs.size() && failure.get()==null; idx=next.getAndIncrement()){
					// Claim coarse per-family work; the main thread joins by dense rank and owns shared writers.
					final Job job=jobs.get(idx);
					final Result result=exportJob(out, job, bgPin);
					assert(results[job.rank]==null) : "Duplicate dense rank from tail-family registry preflight: "+job.rank;
					results[job.rank]=result;
				}
			}catch(Throwable e){failure.compareAndSet(null, e);}
		}
		private final Path out;
		private final ArrayList<Job> jobs;
		private final Result[] results;
		private final String bgPin;
		private final AtomicInteger next;
		private final AtomicReference<Throwable> failure;
	}

	private static final class Job{
		Job(final int rank, final String family, final String retained, final Path directory){
			this.rank=rank; this.family=family; this.retained=retained; this.directory=directory;
		}
		private final int rank;
		private final String family;
		private final String retained;
		private final Path directory;
	}

	private static final class Result{
		Result(final int rank, final long residues, final ByteBuilder consensus, final ByteBuilder family, final ByteBuilder model){
			this.rank=rank; this.residues=residues; this.consensus=consensus; this.family=family; this.model=model;
		}
		private final int rank;
		private final long residues;
		private final ByteBuilder consensus;
		private final ByteBuilder family;
		private final ByteBuilder model;
	}

	private static final String REGISTRY_HEADER="family_id\tinput_copies\tretained_copies\trejected_copies\tstate\tdirectory";
	static final String MODELS_HEADER="rank\tfamily_id\tdirectory\tmembers\tresidues\tbundle_sha80\tconsensus_sha80\tprovenance_sha80\tretained_sha80";
	private static final String[] ARTIFACTS={"refined.mqhb", "consensus.faa", "provenance.tsv", "retained.faa"};
	private HbmTailLibrary(){}
}
