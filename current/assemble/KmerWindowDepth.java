package assemble;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import bloom.KmerCountAbstract;
import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.ListNum;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import stream.Writer;
import stream.WriterFactory;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/**
 * Reports exact read-kmer evidence along selected reference windows.
 * Coordinates describe forward-reference kmer ends, one-based and inclusive.
 * Reference sequence selects queries only; it never contributes counts.
 * @author Fischl
 */
public final class KmerWindowDepth {

	/** Loads one exact table and emits depths without assembling or mutating it. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, KmerWindowDepth.class, false);
		try{new KmerWindowDepth(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Separates reference-window options from the existing native table parser. */
	private KmerWindowDepth(final String[] args){
		assert(args!=null) : "PreParser supplies expanded diagnostic arguments.";
		final ArrayList<String> tableArgs=new ArrayList<String>();
		String ranges=null;
		boolean hollow=false;
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("ref")){reference=b;}
			else if(a.equals("queries")){queries=b;}
			else if(a.equals("querystart")){queryStart=Parse.parseIntKMG(b);}
			else if(a.equals("scaffold")){scaffold=b;}
			else if(a.equals("regions")){ranges=b;}
			else if(a.equals("out")){output=b;}
			else if(a.equals("outseq")){sequenceOutput=b;}
			else if(a.equals("pathseedlen")){pathSeedLength=Parse.parseIntKMG(b);}
			else if(a.equals("pathratio")){pathRatio=Parse.parseIntKMG(b);}
			else if(a.equals("pathsupport")){pathSupport=Parse.parseIntKMG(b);}
			else if(a.equals("pathbothspans")){pathBothSpans=Parse.parseBoolean(b);}
			else if(a.equals("pathrollback")){pathRollback=Parse.parseBoolean(b);}
			else if(a.equals("pathbranchstop")){pathBranchStop=Parse.parseBoolean(b);}
			else if(a.equals("pathselfstop")){pathSelfStop=Parse.parseBoolean(b);}
			else if(a.equals("pathbackstop")){pathBackStop=Parse.parseBoolean(b);}
			else if(a.equals("pathbackratio")){pathBackRatio=Parse.parseIntKMG(b);}
			else if(a.equals("pathbackactualsupport")){pathBackActualSupport=Parse.parseIntKMG(b);}
			else if(a.equals("pathdropstop")){pathDropStop=Parse.parseBoolean(b);}
			else if(a.equals("pathdropratio")){pathDropRatio=Parse.parseIntKMG(b);}
			else if(a.equals("pathdropcontext")){pathDropContext=Parse.parseIntKMG(b);}
			else if(a.equals("pathdropmax")){pathDropMax=Parse.parseIntKMG(b);}
			else if(a.equals("pathrisestop")){pathRiseStop=Parse.parseBoolean(b);}
			else if(a.equals("pathriseratio")){pathRiseRatio=Parse.parseIntKMG(b);}
			else if(a.equals("pathrisemin")){pathRiseMin=Parse.parseIntKMG(b);}
			else if(a.equals("pathlowbranchstop")){pathLowBranchStop=Parse.parseBoolean(b);}
			else if(a.equals("pathlowbranchratio")){pathLowBranchRatio=Parse.parseIntKMG(b);}
			else if(a.equals("pathsolidstop")){pathSolidStop=Parse.parseBoolean(b);}
			else if(a.equals("pathsolidbranchstop")){pathSolidBranchStop=Parse.parseBoolean(b);}
			else if(a.equals("pathnominalstop")){pathNominalStop=Parse.parseBoolean(b);}
			else if(a.equals("pathnominalweak")){pathNominalWeak=Parse.parseBoolean(b);}
			else if(a.equals("pathspacingabstain")){pathSpacingAbstain=Parse.parseBoolean(b);}
			else if(a.equals("pathshortreads")){pathShortReads=b;}
			else if(a.equals("pathshortk")){pathShortK=Parse.parseIntKMG(b);}
			else if(a.equals("pathshortprefilter")){pathShortPrefilter=Parse.parseBoolean(b);}
			else if(a.equals("pathnominalrise")){pathNominalRise=Parse.parseBoolean(b);}
			else if(a.equals("pathrisesource")){pathRiseSource=Parse.parseBoolean(b);}
			else if(a.equals("pathsourcek")){pathSourceK=Parse.parseIntKMG(b);}
			else if(a.equals("pathshortstop")){pathShortStop=Parse.parseBoolean(b);}
			else if(a.equals("pathshortdrop")){pathShortDrop=Parse.parseBoolean(b);}
			else if(a.equals("pathsolidbackstop")){pathSolidBackStop=Parse.parseBoolean(b);}
			else if(a.equals("pathbackfork")){pathBackFork=Parse.parseBoolean(b);}
			else if(a.equals("pathcandidates")){pathCandidates=b;}
			else if(a.equals("pathsolidsupport")){pathSolidSupport=Parse.parseIntKMG(b);}
			else if(a.equals("selfseeds")){selfSeeds=b;}
			else if(a.equals("hollow")){hollow=Parse.parseBoolean(b);}
			else if(a.equals("queryspan")){querySpan=Integer.parseInt(b);}
			else{tableArgs.add(arg);}
		}
		if(output==null || (queries==null && (reference==null || scaffold==null || ranges==null))){
			throw new IllegalArgumentException("Required: in=reads ref=fasta scaffold=name regions=start-end[,start-end] out=tsv k=N");
		}
		if(queries!=null && (reference!=null || scaffold!=null || ranges!=null)){
			throw new IllegalArgumentException("queries= terminal mode is exclusive with ref/scaffold/regions.");
		}
		for(String range : ranges==null ? new String[0] : ranges.split(",", -1)){
			final String[] ends=range.split("-", -1);
			if(ends.length!=2){throw new IllegalArgumentException("Invalid inclusive region: "+range);}
			regions.add(new int[]{Integer.parseInt(ends[0]), Integer.parseInt(ends[1])});
		}
		// Hollow requires paired hashes. Solid retains the original explicit table.
		// The query words below must use the same geometry as the table's loader.
		tableArgs.add("prefilter=f");
		tableArgs.add("prealloc=f");
		tableArgs.add("hashonly="+hollow);
		tableArgs.add("hollow="+hollow);
		tableArgs.add("rcomp=t");
		tableArgs.add("minprob=0");
		tableArgs.add("minprobmain=f");
		tableArgs.add("packed=t");
		if(pathSelfStop){
			//The self ledger counts raw FASTA windows, so load precisely that evidence.
			tableArgs.add("qtrim=f");
			tableArgs.add("maq=0");
			tableArgs.add("ecco=f");
			tableArgs.add("merge=f");
			tableArgs.add("reads=-1");
			tableArgs.add("onepass=f");
		}
		tables=new KmerTableSetU(tableArgs.toArray(new String[0]), 0);
		k=tables.kbig;
		if(querySpan==0){querySpan=k;}
		if(querySpan<1 || (querySpan!=k && (!hollow || querySpan<66 || k<66))){
			throw new IllegalArgumentException("Unequal queryspan and k require hollow=t and both spans>=66 for neighbor context.");
		}
		if(queryStart<0 || (queryStart>0 && (queries==null || queryStart<querySpan || pathSeedLength!=0 ||
				sequenceOutput!=null || pathBothSpans || pathRollback || pathBranchStop || pathSelfStop ||
				pathBackStop || pathDropStop || pathRiseStop || pathLowBranchStop || pathSolidStop || pathSolidBackStop || pathShortStop || pathNominalStop))){
			throw new IllegalArgumentException("querystart must be zero or >=queryspan with queries supplied and all path policies/outseq off.");
		}
		if(pathSeedLength<0 || pathRatio<1 || pathSupport<1){
			throw new IllegalArgumentException("Require pathseedlen>=0, pathratio>=1 and pathsupport>=1.");
		}
		if(pathBackRatio<2 || pathDropRatio<2 || pathDropContext<1 || pathDropMax<0){
			throw new IllegalArgumentException("Require pathbackratio>=2, pathdropratio>=2, pathdropcontext>=1 and pathdropmax>=0.");
		}
		if(pathBackActualSupport<0){
			throw new IllegalArgumentException("Require pathbackactualsupport>=0; zero inherits pathsupport.");
		}
		if(pathBackActualSupport==0){pathBackActualSupport=pathSupport;}
		if(pathRiseRatio<2 || pathRiseMin<1){
			throw new IllegalArgumentException("Require pathriseratio>=2 and pathrisemin>=1.");
		}
		if(pathLowBranchRatio<2){throw new IllegalArgumentException("Require pathlowbranchratio>=2.");}
		if(pathSolidSupport<2){throw new IllegalArgumentException("Require pathsolidsupport>=2; one observation must not establish a solid contradiction.");}
		if(pathSolidBranchStop && (!pathSolidStop || pathRatio<2)){
			throw new IllegalArgumentException("pathsolidbranchstop requires pathsolidstop=t and pathratio>=2.");
		}
		if(pathShortDrop && !pathShortStop){throw new IllegalArgumentException("pathshortdrop requires pathshortstop=t; it replaces that diagnostic's predicate only.");}
		final boolean solidMode=pathSolidStop || pathSolidBackStop;
		if(pathNominalWeak && (!pathNominalStop || pathRatio<2)){
			throw new IllegalArgumentException("pathnominalweak requires pathnominalstop=t and pathratio>=2; only count1 choices become eligible.");
		}
		if(pathSpacingAbstain && (!pathNominalStop || pathNominalWeak || pathNominalRise || k<67 || pathSeedLength<=k)){
			throw new IllegalArgumentException("pathspacingabstain requires pathnominalstop=t, k>=67, pathseedlen>=k+1, pathnominalweak=f and pathnominalrise=f.");
		}
		if(pathNominalRise && !pathNominalStop){
			throw new IllegalArgumentException("pathnominalrise requires pathnominalstop=t; it adds a rise check to that exclusive Hollow scan.");
		}
		if(pathShortReads!=null && (!pathNominalStop || pathNominalWeak || pathNominalRise || pathSourceK!=0)){
			throw new IllegalArgumentException("pathshortreads requires zero-choice pathnominalstop=t, pathnominalweak=f, pathnominalrise=f and pathsourcek=0.");
		}
		if(pathShortPrefilter && (pathShortReads==null || KmerCountAbstract.SKETCH_MODE)){
			throw new IllegalArgumentException("pathshortprefilter requires pathshortreads and complete occurrence counting, not sketch sampling.");
		}
		if((pathShortK!=64 && pathShortK!=96) || (pathShortK!=64 && pathShortReads==null) ||
				(pathShortReads!=null && pathShortK>k)){
			throw new IllegalArgumentException("pathshortk must be64 or96, fit within k, and nondefault96 requires pathshortreads.");
		}
		if(pathShortPrefilter && pathShortK!=64){
			throw new IllegalArgumentException("pathshortprefilter is validated only for pathshortk=64; use unfiltered counts at96.");
		}
		if(pathSourceK<0 || (pathSourceK>0 && (!pathNominalStop || pathSourceK<=64 || pathSourceK>=k))){
			throw new IllegalArgumentException("pathsourcek must be zero, or 64<pathsourcek<k with pathnominalstop=t.");
		}
		if(pathRiseSource && (!pathNominalRise || pathSourceK==0)){
			throw new IllegalArgumentException("pathrisesource requires pathnominalrise=t and a positive pathsourcek.");
		}
		if(pathNominalStop && (!hollow || k<=64 || querySpan!=k || pathSeedLength<k || pathSupport<2 ||
				solidMode || pathShortStop || pathBothSpans || pathRollback || pathBranchStop || pathSelfStop ||
				pathBackStop || pathDropStop || pathRiseStop || pathLowBranchStop)){
			throw new IllegalArgumentException("pathnominalstop requires hollow=t, k>64, queryspan=k, pathseedlen>=k, pathsupport>=2 and every other path policy off.");
		}
		if(pathShortStop && (!hollow || k<=64 || querySpan!=k || pathSeedLength<k || solidMode ||
				pathBothSpans || pathRollback || pathBranchStop || pathSelfStop || pathBackStop ||
				pathDropStop || pathRiseStop || pathLowBranchStop)){
			throw new IllegalArgumentException("pathshortstop requires hollow=t, k>64, queryspan=k, pathseedlen>=k and every other path policy off.");
		}
		if(solidMode && (hollow || querySpan!=k || pathSeedLength<k || pathBothSpans ||
				pathRollback || pathBranchStop || pathSelfStop || pathBackStop || pathDropStop || pathRiseStop || pathLowBranchStop)){
			throw new IllegalArgumentException("Solid path modes require hollow=f, queryspan=k, pathseedlen>=k and all Hollow path policies off.");
		}
		if(pathSolidBackStop && (pathSolidStop || pathCandidates==null || pathSeedLength<=k || k<67 || pathSolidSupport!=2)){
			throw new IllegalArgumentException("pathsolidbackstop requires pathcandidates, k>=67, pathseedlen>k, pathsolidsupport=2 and pathsolidstop=f.");
		}
		if(!pathSolidBackStop && pathCandidates!=null){throw new IllegalArgumentException("pathcandidates requires pathsolidbackstop=t.");}
		if(pathBackFork && !pathSolidBackStop){throw new IllegalArgumentException("pathbackfork requires pathsolidbackstop=t.");}
		if((pathBackStop || pathDropStop || pathRiseStop || pathLowBranchStop) && pathSeedLength<=k){
			throw new IllegalArgumentException("Nominal path checks require pathseedlen>=k+1, including the actual predecessor base.");
		}
		if(pathBothSpans && (pathSeedLength==0 || querySpan!=k || k<67)){
			throw new IllegalArgumentException("pathbothspans requires positive pathseedlen, k>=67 and queryspan omitted or equal to k.");
		}
		if(pathRollback && pathSeedLength==0){
			throw new IllegalArgumentException("pathrollback requires positive pathseedlen; it never edits the original seed.");
		}
		if(pathBranchStop && (pathSeedLength==0 || pathRollback)){
			throw new IllegalArgumentException("pathbranchstop requires positive pathseedlen and pathrollback=f; these are separate diagnostic policies.");
		}
		pathSpans=pathBothSpans ? new int[]{k+1, k-1} : new int[]{querySpan};
		if(pathSeedLength>0 && (queries==null || sequenceOutput==null || pathSeedLength<pathSpans[0])){
			throw new IllegalArgumentException("Path veto requires queries=extensions, outseq=guarded.fa and a seed at least as long as every query span.");
		}
		if(pathSeedLength>0 && !solidMode && !pathShortStop && !pathNominalStop && (!hollow || querySpan<66 || k<66 ||
				(!pathBothSpans && Math.abs(k-querySpan)!=1))){
			throw new IllegalArgumentException("Path veto requires queries=extensions outseq=guarded.fa, hollow=t, spans>=66 and a seed at least as long as every query span. Use k=queryspan+/-1 or pathbothspans=t.");
		}
		if(pathSeedLength==0 && sequenceOutput!=null){
			throw new IllegalArgumentException("outseq is only valid with positive pathseedlen.");
		}
		final ArrayList<String> inputs=new ArrayList<String>(tables.in1);
		inputs.addAll(tables.in2);
		inputs.addAll(tables.extra);
		if(pathSelfStop){
			if(pathSeedLength<k || selfSeeds==null || !tables.in2.isEmpty()){
				throw new IllegalArgumentException("pathselfstop requires selfseeds=original.fa, pathseedlen>=k and no paired in2.");
			}
			if(!FileFormat.testInput(selfSeeds, FileFormat.FASTA, null, true, true).fasta() ||
					!FileFormat.testInput(queries, FileFormat.FASTA, null, true, true).fasta()){
				throw new IllegalArgumentException("pathselfstop requires FASTA originals and saved extensions.");
			}
			int copies=0;
			for(String input : inputs){if(sameFile(input, selfSeeds)){copies++;}}
			if(copies!=1){throw new IllegalArgumentException("selfseeds must occur exactly once in in/extra; found "+copies);}
		}else if(selfSeeds!=null){throw new IllegalArgumentException("selfseeds requires pathselfstop=t.");}
		final String queryFile=queries==null ? reference : queries;
		if(pathShortReads!=null){
			int copies=0;
			for(String input : inputs){if(sameFile(input, pathShortReads)){copies++;}}
			if(copies!=1 || !tables.in2.isEmpty() || sameFile(pathShortReads, queryFile) ||
					!FileFormat.testInput(pathShortReads, FileFormat.FASTQ, null, true, true).fastq()){
				throw new IllegalArgumentException("pathshortreads requires no in2 and one original FASTQ counted exactly once in in/extra, never the query or seed FASTA.");
			}
		}
		if(pathCandidates!=null){
			if(sameFile(pathCandidates, queryFile)){throw new IllegalArgumentException("Candidates and saved sequences must be separate inputs.");}
			for(String input : inputs){
				if(sameFile(input, pathCandidates)){throw new IllegalArgumentException("Candidate decisions must never contribute kmer counts.");}
			}
			inputs.add(pathCandidates);//Protect the report from output aliases without counting it.
		}
		TadpoleGraph.checkPaths(queryFile, inputs, output, sequenceOutput);
		if(!Tools.testInputFiles(false, true, queryFile, pathCandidates) ||
				!Tools.testOutputFiles(false, false, false, output, sequenceOutput)){
			throw new IllegalArgumentException("Missing reference or existing/unwritable output.");
		}
		if(pathShortReads!=null){
			//Do not inherit seed inputs or read-sampling/trimming filters: this table
			//counts only the explicitly named complete original FASTQ, as in the profile audit.
			final String[] shortArgs={"in="+pathShortReads, "k="+pathShortK, "hollow=f", "hashonly=f",
					"packed=t", "rcomp=t", "prefilter="+(pathShortPrefilter ? "1" : "f"),
					"prepasses=1", "prehashes=2", "prealloc=f", "minprob=0", "minprobmain=f",
					"qtrim=f", "maq=0", "ecco=f", "merge=f", "reads=-1", "onepass=f",
					"interleaved=f", "t="+Shared.threads()};
			shortTables=new KmerTableSetU(shortArgs, 0);
		}else if(pathShortStop || pathBackFork){
			//Both constructors use the same packed/canonical/probability settings.
			//Only the second table's geometry changes; neither table sees queries.
			tableArgs.add("k=64");
			tableArgs.add("hollow=f");
			tableArgs.add("hashonly=f");
			shortTables=new KmerTableSetU(tableArgs.toArray(new String[0]), 0);
		}else{shortTables=null;}
		if(pathSourceK>0){
			//The shorter source span counts exactly the same originals, never the saved queries.
			//This mode is exclusive with all short/Solid policies, so no third table is loaded.
			tableArgs.add("k="+pathSourceK);
			tableArgs.add("hollow=t");
			tableArgs.add("hashonly=t");
			sourceTables=new KmerTableSetU(tableArgs.toArray(new String[0]), 0);
			assert(sourceTables.kbig==pathSourceK && sourceTables.hollow) :
					"Source-choice queries must match the explicitly requested packed Hollow geometry.";
		}else{sourceTables=null;}
	}

	/** Resolves aliases before checking that the original seed input is counted once. */
	private static boolean sameFile(final String a, final String b){
		assert(a!=null && b!=null) : "The evidence inclusion check requires two real file names.";
		try{return new File(a).getCanonicalFile().equals(new File(b).getCanonicalFile());}
		catch(IOException e){throw new IllegalArgumentException("Could not resolve seed evidence paths.", e);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Validates every requested base before loading reads or opening output. */
	private byte[] referenceBases(){
		byte[] bases=null;
		for(Contig contig : TadpoleGraph.readContigs(reference)){
			final String name=contig.name;
			int length=0;
			while(length<name.length() && !Character.isWhitespace(name.charAt(length))){length++;}
			if(!scaffold.equals(name.substring(0, length))){continue;}
			if(bases!=null){throw new IllegalArgumentException("Duplicate scaffold identifier: "+scaffold);}
			bases=contig.bases;
		}
		if(bases==null){throw new IllegalArgumentException("Scaffold not found: "+scaffold);}
		for(int[] region : regions){
			if(region[0]<querySpan || region[1]<region[0] || region[1]>bases.length){
				throw new IllegalArgumentException("Invalid kmer-end region "+region[0]+"-"+region[1]+" for queryspan="+querySpan+", length="+bases.length);
			}
			for(int pos=region[0]-querySpan; pos<region[1]; pos++){
				if(AminoAcid.baseToNumber[bases[pos]]<0){
					throw new IllegalArgumentException("Ambiguous reference base at "+scaffold+":"+(pos+1));
				}
			}
		}
		return bases;
	}

	/** Owns writer and table cleanup even when parsing or output fails. */
	private void process(){
		if(pathShortStop){processShortPaths(); return;}
		if(queryStart>0){processQueryWindows(); return;}
		if(pathSeedLength>0){processPaths(); return;}
		if(queries!=null){processTerminals(); return;}
		final byte[] bases=referenceBases();
		ByteStreamWriter writer=null;
		boolean complete=false, writeError=false;
		long rows=0;
		try{
			tables.process(new Timer());
			writer=new ByteStreamWriter(output, false, false, false);
			writer.start();
			writer.print("region\tscaffold\tk"+(querySpan==k ? "" : "\tquery_span")+
					"\tstart1\tend1\tword\tdepth\tleft_A\tleft_C\tleft_G\tleft_T\tright_A\tright_C\tright_G\tright_T\tref_previous\tref_next\n");
			final Kmer word=new Kmer(k, tables.hollow), neighbor=new Kmer(k, tables.hollow);
			final ByteBuilder row=new ByteBuilder(k+200);
			for(int index=0; index<regions.size(); index++){
				final int[] region=regions.get(index);
				word.clear();
				if(querySpan==k){
					for(int pos=region[0]-k; pos<region[0]-1; pos++){
						word.addRightNumeric(AminoAcid.baseToNumber[bases[pos]]);
					}
				}
				for(int end=region[0]; end<=region[1]; end++){
					if(querySpan==k){word.addRightNumeric(AminoAcid.baseToNumber[bases[end-1]]);}
					else{setHollowQuery(word, bases, end-querySpan, querySpan);}
					assert(word.len>=k) : "Each output row must contain a complete query word.";
					row.clear().append(index+1).tab().append(scaffold).tab().append(k);
					if(querySpan!=k){row.tab().append(querySpan);}
					row.tab().append(end-querySpan+1).tab().append(end).tab().append(bases, end-querySpan, querySpan);
					row.tab().append(depth(word));
					for(int base=0; base<4; base++){
						neighbor.setFrom(word); neighbor.addLeftNumeric(base);
						row.tab().append(depth(neighbor));
					}
					for(int base=0; base<4; base++){
						neighbor.setFrom(word); neighbor.addRightNumeric(base);
						row.tab().append(depth(neighbor));
					}
					row.tab().append(end>querySpan ? (char)bases[end-querySpan-1] : '.');
					row.tab().append(end<bases.length ? (char)bases[end] : '.').nl();
					writer.print(row);
					rows++;
				}
			}
			complete=true;
		}finally{
			if(writer!=null){writeError=writer.poisonAndWait();}
			tables.clear();
			if(writeError && !complete){System.err.println("Window output also failed during diagnostic cleanup.");}
		}
		if(writeError){throw new RuntimeException("Window output failed: "+output);}
		System.err.println("KMER_WINDOW_DEPTH_PASS rows="+rows+" k="+k);
	}

	/** Reports every saved next-base context from a fixed boundary, loading evidence once.
	 * Coordinates and counts use the saved sequence orientation, not canonical orientation.
	 * The last row precedes the final saved base; a record ending at the boundary has no rows.
	 * This diagnostic neither changes bases nor decides where extension should stop. */
	private void processQueryWindows(){
		assert(queries!=null && queryStart>=querySpan && pathSeedLength==0) :
				"Constructor validation makes every initial context complete and excludes sequence-changing path modes.";
		final Streamer stream=StreamerFactory.makeStreamer(
				FileFormat.testInput(queries, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 0);
		ByteStreamWriter writer=null;
		boolean started=false, error=false;
		long reads=0, rows=0;
		try{
			tables.process(new Timer());
			writer=new ByteStreamWriter(output, false, false, false);
			writer.start();
			writer.print("id\tk\tquery_span\tstart1\tend1\tword\tdepth\tleft_A\tleft_C\tleft_G\tleft_T"+
					"\tright_A\tright_C\tright_G\tright_T\tquery_previous\tquery_next\n");
			final Kmer word=new Kmer(k, tables.hollow), neighbor=new Kmer(k, tables.hollow);
			final ByteBuilder row=new ByteBuilder(querySpan+200);
			stream.start(); started=true;
			for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
				for(Read read : batch){
					if(read.mate!=null || read.length()<queryStart || read.id==null || read.id.indexOf('\t')>=0){
						throw new IllegalArgumentException("Window queries require unpaired records at least querystart long with tab-free IDs: "+read.id);
					}
					final byte[] bases=read.bases;
					for(int pos=queryStart-querySpan; pos<read.length(); pos++){
						if(AminoAcid.baseToNumber[bases[pos]]<0){throw new IllegalArgumentException("Ambiguous window query: "+read.id+" at "+(pos+1));}
					}
					word.clear();
					if(querySpan==k){
						for(int pos=queryStart-k; pos<queryStart-1; pos++){word.addRightNumeric(AminoAcid.baseToNumber[bases[pos]]);}
					}
					for(int end=queryStart; end<read.length(); end++){
						if(querySpan==k){word.addRightNumeric(AminoAcid.baseToNumber[bases[end-1]]);}
						else{setHollowQuery(word, bases, end-querySpan, querySpan);}
						row.clear().append(read.id).tab().append(k).tab().append(querySpan)
								.tab().append(end-querySpan+1).tab().append(end).tab().append(bases, end-querySpan, querySpan)
								.tab().append(depth(word));
						for(int base=0; base<4; base++){
							neighbor.setFrom(word); neighbor.addLeftNumeric(base);
							row.tab().append(depth(neighbor));
						}
						for(int base=0; base<4; base++){
							neighbor.setFrom(word); neighbor.addRightNumeric(base);
							row.tab().append(depth(neighbor));
						}
						row.tab().append(end>querySpan ? (char)bases[end-querySpan-1] : '.');
						row.tab().append((char)bases[end]).nl();
						writer.print(row);
						rows++;
					}
					reads++;
				}
			}
		}finally{
			if(writer!=null){error|=writer.poisonAndWait();}
			if(started){error|=ReadWrite.closeStream(stream);}
			tables.clear();
		}
		if(error){throw new RuntimeException("Batch window query I/O failed.");}
		System.err.println("KMER_QUERY_WINDOWS_PASS reads="+reads+" rows="+rows+" k="+k+" queryspan="+querySpan);
	}

	/** Streams one terminal window per single-end seed; queries never add counts. */
	private void processTerminals(){
		assert(queries!=null && querySpan>=64) : "Terminal reports include two32bp anchor fields for downstream identity checks.";
		if(querySpan<64){throw new IllegalArgumentException("Terminal queries require queryspan>=64.");}
		final FileFormat ff=FileFormat.testInput(queries, FileFormat.FASTA, null, true, true);
		final Streamer stream=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 1);
		ByteStreamWriter writer=null;
		boolean started=false, error=false;
		long rows=0;
		try{
			tables.process(new Timer());
			writer=new ByteStreamWriter(output, false, false, false);
			writer.start();
			writer.print("id\tlength\tk\tquery_span\tleft_anchor\tright_anchor\tdepth\n");
			final Kmer word=new Kmer(k, tables.hollow);
			final ByteBuilder row=new ByteBuilder();
			stream.start(); started=true;
			for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
				for(Read read : batch){
					if(read.mate!=null || read.length()<querySpan || read.id.indexOf('\t')>=0){
						throw new IllegalArgumentException("Terminal query requires a long-enough unpaired seed and a tab-free ID: "+read.id);
					}
					final int start=read.length()-querySpan;
					for(int i=start; i<read.length(); i++){
						if(AminoAcid.baseToNumber[read.bases[i]]<0){throw new IllegalArgumentException("Ambiguous terminal query: "+read.id);}
					}
					if(querySpan!=k){setHollowQuery(word, read.bases, start, querySpan);}
					else{
						word.clear();
						for(int i=start; i<read.length(); i++){word.addRight(read.bases[i]);}
					}
					row.clear().append(read.id).tab().append(read.length()).tab().append(k).tab().append(querySpan)
							.tab().append(read.bases, start, 32).tab().append(read.bases, read.length()-32, 32)
							.tab().append(depth(word)).nl();
					writer.print(row);
					rows++;
				}
			}
		}finally{
			if(writer!=null){error|=writer.poisonAndWait();}
			if(started){error|=ReadWrite.closeStream(stream);}
			tables.clear();
		}
		if(error){throw new RuntimeException("Terminal query I/O failed.");}
		System.err.println("KMER_TERMINAL_DEPTH_PASS rows="+rows+" k="+k+" queryspan="+querySpan);
	}

	/** Tests Hollow choices against independently counted Solid64 alternatives.
	 * Both tables count the same original evidence once. Queries are streamed unchanged
	 * except for a rejected added suffix; no replacement base or reference is used. */
	private void processShortPaths(){
		assert(pathShortStop && tables.hollow && shortTables!=null && !shortTables.hollow && shortTables.kbig==64) :
				"The exclusive short-key guard requires distinct Hollow and explicit Solid64 evidence tables.";
		final Streamer stream=StreamerFactory.makeStreamer(
				FileFormat.testInput(queries, FileFormat.FASTA, null, true, true), null, true, -1, false, true, 0);
		Writer sequences=null;
		ByteStreamWriter report=null;
		boolean streamStarted=false, writerStarted=false, error=false;
		long reads=0, trimmed=0, removed=0;
		try{
			tables.process(new Timer());
			shortTables.process(new Timer());
			report=new ByteStreamWriter(output, false, false, false);
			report.start();
			report.print("id\tseed_length\toriginal_length\tretained_length\tk\tshort_k\thollow_chosen\tsolid_chosen\tsolid_other\tother_base\tstop_reason"+
					(pathShortDrop ? "\thollow_context\tsolid_context" : "")+"\n");
			sequences=WriterFactory.makeWriter(
					FileFormat.testOutput(sequenceOutput, FileFormat.FASTA, null, true, false, false, true), null, 1, null, false);
			sequences.start(); writerStarted=true;
			final Kmer word=new Kmer(k, true), shortWord=new Kmer(64, false), neighbor=new Kmer(64, false);
			final int[] evidence=new int[pathShortDrop ? 6 : 4];
			final ByteBuilder row=new ByteBuilder();
			stream.start(); streamStarted=true;
			for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
				for(Read read : batch){
					validatePath(read);
					final int original=read.length();
					final int retained=pathShortDrop ? shortDropLimit(read, word, shortWord, neighbor, evidence) :
							shortLimit(read, word, shortWord, neighbor, evidence);
					assert(retained>=pathSeedLength && retained<=original) : "A short-key veto must preserve the original seed and cannot grow a saved path.";
					row.clear().append(read.id).tab().append(pathSeedLength).tab().append(original).tab().append(retained)
							.tab().append(k).tab().append(64).tab().append(evidence[0]).tab().append(evidence[1])
							.tab().append(evidence[2]).tab().append(evidence[3]<0 ? '.' : (char)AminoAcid.numberToBase[evidence[3]])
							.tab().append(retained<original ? (pathShortDrop ? "hollow_drop_solid_contradiction" : "weak_hollow_short_dominance") : "none");
					if(pathShortDrop){row.tab().append(evidence[4]).tab().append(evidence[5]);}
					row.nl();
					report.print(row);
					if(retained<original){
						trimmed++; removed+=original-retained;
						read.bases=Arrays.copyOf(read.bases, retained);
						if(read.quality!=null){read.quality=Arrays.copyOf(read.quality, retained);}
					}
					reads++;
				}
				sequences.addReads(batch);
			}
		}finally{
			if(report!=null){error|=report.poisonAndWait();}
			if(writerStarted){error|=sequences.poisonAndWait();}
			if(streamStarted){error|=ReadWrite.closeStream(stream);}
			tables.clear(); shortTables.clear();
		}
		if(error){throw new RuntimeException("Short-key path sequence/report I/O failed.");}
		System.err.println("KMER_SHORT_PATH_PASS reads="+reads+" trimmed="+trimmed+" removed_bases="+removed+" k="+k+" shortk=64");
	}

	/** Advances both observed geometries; only a weak Hollow next key needs short branching.
	 * Scratch is restored from the short context before each candidate. The returned
	 * boundary precedes the conflicting saved base, never an earlier seed position. */
	private int shortLimit(final Read read, final Kmer word, final Kmer shortWord, final Kmer neighbor, final int[] evidence){
		assert(read.length()>=pathSeedLength && pathSeedLength>=k && evidence.length==4) :
				"Validated paths contain both complete seed contexts and four diagnostic evidence slots.";
		Arrays.fill(evidence, -1);
		word.clear(); shortWord.clear();
		for(int i=pathSeedLength-k; i<pathSeedLength; i++){word.addRight(read.bases[i]);}
		for(int i=pathSeedLength-64; i<pathSeedLength; i++){shortWord.addRight(read.bases[i]);}
		for(int end=pathSeedLength; end<read.length(); end++){
			final int chosenBase=AminoAcid.baseToNumber[read.bases[end]];
			word.addRightNumeric(chosenBase);
			final int hollowChosen=depth(word);
			if(hollowChosen<3){
				neighbor.setFrom(shortWord); neighbor.addRightNumeric(chosenBase);
				final int solidChosen=Math.max(0, shortTables.getCount(neighbor));
				if(solidChosen>=3){
					int other=0, otherBase=-1;
					for(int base=0; base<4; base++){
						if(base==chosenBase){continue;}
						neighbor.setFrom(shortWord); neighbor.addRightNumeric(base);
						final int count=Math.max(0, shortTables.getCount(neighbor));
						if(count>other){other=count; otherBase=base;}
					}
					if(rejectShortStep(hollowChosen, solidChosen, other)){
						evidence[0]=hollowChosen; evidence[1]=solidChosen; evidence[2]=other; evidence[3]=otherBase;
						return end;
					}
				}
			}
			shortWord.addRightNumeric(chosenBase);
		}
		return read.length();
	}

	/** Fixed experimental conjunction; neither weak long evidence nor a short branch suffices alone. */
	public static boolean rejectShortStep(final int hollowChosen, final int solidChosen, final int solidOther){
		assert(hollowChosen>=0 && solidChosen>=0 && solidOther>=0) :
				"Both exact tables normalize absent keys to zero before testing short-key dominance.";
		return hollowChosen<3 && solidChosen>=3 && solidOther>=2L*solidChosen;
	}

	/** Requires a supported Hollow collapse and an independently supported short contradiction.
	 * Counts roll along the complete saved suffix, including positions after an earlier
	 * short contradiction that did not collapse. The original seed always remains intact. */
	private int shortDropLimit(final Read read, final Kmer word, final Kmer shortWord,
			final Kmer neighbor, final int[] evidence){
		assert(pathShortDrop && tables.hollow && shortTables!=null && evidence.length==6 &&
				pathSeedLength>=k && read.length()>=pathSeedLength) :
				"The two-table diagnostic requires complete seed contexts and six evidence slots; only added bases may be removed.";
		Arrays.fill(evidence, -1);
		word.clear(); shortWord.clear();
		for(int i=pathSeedLength-k; i<pathSeedLength; i++){word.addRight(read.bases[i]);}
		for(int i=pathSeedLength-64; i<pathSeedLength; i++){shortWord.addRight(read.bases[i]);}
		int hollowContext=depth(word);
		for(int end=pathSeedLength; end<read.length(); end++){
			final int chosenBase=AminoAcid.baseToNumber[read.bases[end]];
			word.addRightNumeric(chosenBase);
			final int hollowChosen=depth(word);
			if(hollowContext>=32 && hollowContext>=4L*hollowChosen){
				final int solidContext=Math.max(0, shortTables.getCount(shortWord));
				neighbor.setFrom(shortWord); neighbor.addRightNumeric(chosenBase);
				final int solidChosen=Math.max(0, shortTables.getCount(neighbor));
				if(solidContext>=2 && solidChosen==0){
					int other=0, otherBase=-1;
					for(int base=0; base<4; base++){
						if(base==chosenBase){continue;}
						neighbor.setFrom(shortWord); neighbor.addRightNumeric(base);
						final int count=Math.max(0, shortTables.getCount(neighbor));
						if(count>other){other=count; otherBase=base;}
					}
					if(rejectShortDrop(hollowContext, hollowChosen, solidContext, solidChosen, other)){
						evidence[0]=hollowChosen; evidence[1]=solidChosen; evidence[2]=other;
						evidence[3]=otherBase; evidence[4]=hollowContext; evidence[5]=solidContext;
						return end;
					}
				}
			}
			shortWord.addRightNumeric(chosenBase);
			hollowContext=hollowChosen;
		}
		return read.length();
	}

	/** Fixed experimental thresholds: a fourfold drop from depth32 and a Solid64 contradiction.
	 * A flat Hollow path, missing short coverage, or any observed short chosen key abstains. */
	public static boolean rejectShortDrop(final int hollowContext, final int hollowChosen,
			final int solidContext, final int solidChosen, final int solidOther){
		assert(hollowContext>=0 && hollowChosen>=0) :
				"Absent Hollow keys are normalized to zero before the overflow-safe transition comparison.";
		return hollowContext>=32 && hollowContext>=4L*hollowChosen &&
				rejectSolidStep(solidContext, solidChosen, solidOther, 2);
	}

	/** Streams saved extensions through one table, preserving every seed prefix.
	 * Both-span mode tests two observed spacings without reloading counts or query records.
	 * Only the earliest supported next-base conflict is reported, avoiding per-base output.
	 * Solid and nominal Hollow contradiction modes never require blanket coverage.
	 * These experimental vetoes never choose a replacement base or add query counts. */
	private void processPaths(){
		assert(pathSeedLength>=pathSpans[0] && queries!=null && sequenceOutput!=null) :
				"Constructor validation protects the fixed seed boundary and separate query/output streams.";
		final Streamer stream=StreamerFactory.makeStreamer(
				FileFormat.testInput(queries, FileFormat.FASTA, null, true, true), null, true, -1, false, true, pathSelfStop || pathSolidBackStop ? 0 : 1);
		final SeedCursor seeds=pathSelfStop ? new SeedCursor(selfSeeds) : null;
		final HollowSeedCounts self=pathSelfStop ? new HollowSeedCounts(k) : null;
		CandidateCursor candidates=null;
		Writer sequences=null;
		ByteStreamWriter report=null;
		boolean streamStarted=false, writerStarted=false, error=false;
		long reads=0, trimmed=0, removed=0;
		try{
			if(pathSolidBackStop){candidates=new CandidateCursor(pathCandidates);}
			tables.process(new Timer());
			if(pathBackFork || pathShortReads!=null){shortTables.process(new Timer());}
			if(pathShortPrefilter){
				System.err.println("Short-support singleton filtering requested: short_chosen_filtered=0 may mean original depth0 or1; counts>=2 remain exact.");
			}
			if(sourceTables!=null){sourceTables.process(new Timer());}
			report=new ByteStreamWriter(output, false, false, false);
			report.start();
			report.print("id\tseed_length\toriginal_length\tretained_length\tk\tquery_span\tcontext_depth\tchosen_count\tother_count\tother_base"+
					(pathRollback ? "\tconflict_end" : "")+(pathSelfStop ? "\tself_count\texternal_count" : "")+
					(pathBackStop || pathDropStop || pathRiseStop || pathLowBranchStop || pathSolidStop || pathSolidBackStop || pathNominalStop ? "\tstop_reason" : "")+
					(pathSolidBackStop ? "\tsolid_seed_depth" : "")+
					(pathBackFork ? "\tshort_context\tshort_chosen\tshort_other\tshort_other_base\tfork_end" : "")+
					(sourceTables!=null ? "\tsource_k\tsource_chosen" : "")+
					(pathShortReads!=null ? "\tshort_low_windows\t"+(pathShortPrefilter ? "short_chosen_filtered" : "short_chosen") : "")+
					(pathShortReads!=null && pathShortK!=64 ? "\tshort_k" : "")+"\n");
			sequences=WriterFactory.makeWriter(
					FileFormat.testOutput(sequenceOutput, FileFormat.FASTA, null, true, false, false, true), null, 1, null, false);
			sequences.start(); writerStarted=true;
			final Kmer word=new Kmer(k, tables.hollow), neighbor=new Kmer(k, tables.hollow);
			final int[] evidence=new int[7], ambiguity=new int[pathSpans.length];
			final Kmer sourceWord=sourceTables==null ? null : new Kmer(pathSourceK, true);
			final Kmer spacingWord=pathSpacingAbstain ? new Kmer(k, true) : null;
			final ShortKmerPathSupport shortSupport=pathShortReads==null ? null : new ShortKmerPathSupport(shortTables, k);
			final int[] sourceEvidence=sourceTables==null ? null : new int[1];
			final Kmer shortWord=pathBackFork ? new Kmer(64, false) : null;
			final Kmer shortNeighbor=pathBackFork ? new Kmer(64, false) : null;
			final int[] fork=pathBackFork ? new int[5] : null;
			final ByteBuilder row=new ByteBuilder();
			stream.start(); streamStarted=true;
			if(seeds!=null){seeds.start();}
			for(ListNum<Read> batch=stream.nextList(); batch!=null; batch=stream.nextList()){
				for(Read read : batch){
					validatePath(read);
					final int originalLength=read.length();
					int retained;
					if(pathSolidBackStop){retained=solidBackwardLimit(read, word, evidence, candidates, shortWord, shortNeighbor, fork);}
					else if(pathSolidStop || pathNominalStop){retained=contradictionLimit(read, word, neighbor, evidence, sourceWord, sourceEvidence, spacingWord, shortSupport);}
					else{retained=pathLimit(read, word, neighbor, evidence, ambiguity);}
					if(pathBackStop || pathDropStop || pathRiseStop || pathLowBranchStop){retained=nominalLimit(read, retained, word, neighbor, evidence);}
					int selfCount=-1;
					if(seeds!=null){
						final Read seed=seeds.next();
						validateSeed(seed, read);
						self.reset(seed.bases, seed.length());
						word.clear();
						for(int i=pathSeedLength-k; i<pathSeedLength; i++){word.addRight(read.bases[i]);}
						for(int end=pathSeedLength; end<retained; end++){
							neighbor.setFrom(word); neighbor.addRight(read.bases[end]);
							final int global=depth(neighbor), local=self.count(neighbor);
							if(global<local){throw new IllegalStateException("Seed contribution exceeds global count for "+read.id+": "+local+">"+global);}
							if(local>0 && global==local){
								Arrays.fill(evidence, -1);
								evidence[0]=depth(word); evidence[1]=global; evidence[4]=k; evidence[5]=end;
								selfCount=local; retained=end;
								evidence[6]=SELF_ONLY;
								break;
							}
							word.setFrom(neighbor);
						}
					}
					assert(retained>=pathSeedLength && retained<=originalLength) :
							"The path veto may discard added suffix bases only, never the seed or a future base.";
					row.clear().append(read.id).tab().append(pathSeedLength).tab().append(originalLength)
							.tab().append(retained).tab().append(k).tab().append(evidence[4])
							.tab().append(evidence[0]).tab().append(evidence[1]).tab().append(evidence[2])
							.tab().append(evidence[3]<0 ? '.' : (char)AminoAcid.numberToBase[evidence[3]]);
					if(pathRollback){row.tab().append(evidence[5]);}
					if(pathSelfStop){row.tab().append(selfCount).tab().append(selfCount<0 ? -1 : 0);}
					if(pathBackStop || pathDropStop || pathRiseStop || pathLowBranchStop || pathSolidStop || pathSolidBackStop || pathNominalStop){row.tab().append(evidence[6]<0 ? "none" : STOP_REASONS[evidence[6]]);}
					if(candidates!=null){row.tab().append(candidates.solidDepth);}
					if(pathBackFork){
						row.tab().append(fork[0]).tab().append(fork[1]).tab().append(fork[2])
								.tab().append(fork[3]<0 ? '.' : (char)AminoAcid.numberToBase[fork[3]]).tab().append(fork[4]);
					}
					if(sourceTables!=null){row.tab().append(pathSourceK).tab().append(sourceEvidence[0]);}
					if(shortSupport!=null){
						row.tab().append(retained<originalLength ? shortSupport.lowWindowCount() : -1)
							.tab().append(retained<originalLength ? shortSupport.chosenCount() : -1);
						if(pathShortK!=64){row.tab().append(pathShortK);}
					}
					row.nl();
					report.print(row);
					if(retained<originalLength){
						trimmed++; removed+=originalLength-retained;
						read.bases=Arrays.copyOf(read.bases, retained);
						if(read.quality!=null){read.quality=Arrays.copyOf(read.quality, retained);}
					}
					reads++;
				}
				sequences.addReads(batch);
			}
			if(seeds!=null && seeds.next()!=null){throw new IllegalArgumentException("More original seeds than saved extensions.");}
			if(candidates!=null){candidates.requireEnd();}
		}finally{
			if(report!=null){error|=report.poisonAndWait();}
			if(writerStarted){error|=sequences.poisonAndWait();}
			if(streamStarted){error|=ReadWrite.closeStream(stream);}
			if(seeds!=null){error|=seeds.close();}
			if(candidates!=null){error|=candidates.close();}
			tables.clear();
			if(shortTables!=null){shortTables.clear();}
			if(sourceTables!=null){sourceTables.clear();}
		}
		if(error){throw new RuntimeException("Path-veto sequence/report I/O failed.");}
		System.err.println("KMER_PATH_VETO_PASS reads="+reads+" trimmed="+trimmed+" removed_bases="+removed+" k="+k+
				" queryspan="+(pathBothSpans ? (k+1)+","+(k-1) : Integer.toString(querySpan)));
	}

	/** Requires a raw, unpaired FASTA seed matching this saved output exactly. */
	private void validateSeed(final Read seed, final Read extended){
		assert(extended.length()>=pathSeedLength) : "Path validation must precede the exact original-prefix check.";
		if(seed==null || seed.mate!=null || seed.quality!=null || seed.length()!=pathSeedLength || !extended.id.equals(seed.id)){
			throw new IllegalArgumentException("Missing, paired, quality-bearing or reordered original seed for "+extended.id);
		}
		for(int i=0; i<pathSeedLength; i++){
			if(seed.bases[i]!=extended.bases[i]){throw new IllegalArgumentException("Saved extension changed its original seed: "+extended.id);}
		}
	}

	/** Flattens existing streamer batches without loading the seed library into memory. */
	private static final class SeedCursor {
		/** Host-driven FASTA parsing cannot strand a producer on a full record queue after a mismatch. */
		SeedCursor(final String path){
			stream=StreamerFactory.makeStreamer(FileFormat.testInput(path, FileFormat.FASTA, null, true, true),
					null, true, -1, false, true, 0);
		}
		void start(){stream.start(); started=true;}
		/** Returns the next original record, retaining terminal state after the first EOF. */
		Read next(){
			assert(started) : "Original seeds must be streamed in the same order as saved extensions.";
			if(finished){return null;}
			while(batch==null || index>=batch.list.size()){
				batch=stream.nextList(); index=0;
				if(batch==null){finished=true; return null;}
			}
			return batch.list.get(index++);
		}
		boolean close(){return started && ReadWrite.closeStream(stream);}
		private final Streamer stream;
		private ListNum<Read> batch;
		private int index;
		private boolean started, finished;
	}

	/** Streams an existing Hollow actual-floor1 decision table, never count evidence. */
	private static final class CandidateCursor {
		/** Validates the native self/backward report schema before any output is opened. */
		CandidateCursor(final String path){
			file=ByteFile.makeByteFile(path, true);
			boolean valid=false;
			try{
				final byte[] header=file.nextLine();
				valid=header!=null && new String(header, java.nio.charset.StandardCharsets.US_ASCII).equals(CANDIDATE_SCHEMA);
				if(!valid){throw new IllegalArgumentException("Unexpected Hollow candidate schema: "+path);}
			}finally{if(!valid){file.close();}}
		}

		/** Joins by exact ID/order and validates raw-run bounds before selecting seed-start stops.
		 * Original raw extensions may be longer than the already-guarded saved queries. */
		boolean consume(final Read read, final int seedLength, final int tableK, final int[] evidence){
			assert(evidence.length==7 && read.length()>=seedLength) :
					"Path validation and the seven-slot report layout precede candidate consumption.";
			solidDepth=-1;
			final byte[] line=file.nextLine();
			if(line==null){throw new IllegalArgumentException("Hollow candidates ended before query: "+read.id);}
			lp.set(line);
			if(lp.terms()!=13 || !lp.termEquals(read.id, 0)){
				throw new IllegalArgumentException("Hollow candidate schema/ID/order mismatch: "+read.id);
			}
			final int original=lp.parseInt(2), retained=lp.parseInt(3), span=lp.parseInt(5);
			if(lp.parseInt(1)!=seedLength || lp.parseInt(4)!=tableK || original<read.length() ||
					retained<seedLength || retained>original || (span!=-1 && Math.abs((long)span-tableK)>1)){
				throw new IllegalArgumentException("Hollow candidate seed/K/span/length mismatch: "+read.id);
			}
			final int context=lp.parseInt(6), actual=lp.parseInt(7), other=lp.parseInt(8);
			final int self=lp.parseInt(10), external=lp.parseInt(11);
			int otherBase=-1;
			for(int base=0; base<4; base++){
				if(lp.termEquals(AminoAcid.numberToBase[base], 9)){otherBase=base;}
			}
			int reason=-1;
			for(int i=ALTERNATE_SPAN; i<=LOW_FORWARD_BRANCH; i++){
				if(lp.termEquals(STOP_REASONS[i], 12)){reason=i;}
			}
			if(context< -1 || actual< -1 || other< -1 || self< -1 || external< -1 ||
					(otherBase<0 && !lp.termEquals('.', 9)) ||
					(reason<0 && (!lp.termEquals("none", 12) || retained!=original)) ||
					(reason>=0 && (retained>=original || span<1))){
				throw new IllegalArgumentException("Malformed Hollow candidate evidence or stop reason: "+read.id);
			}
			if(reason!=BACKWARD_DOMINANCE){return false;}
			if(span!=tableK || context<0 || actual<0 || other<0 || otherBase<0){
				throw new IllegalArgumentException("Backward candidate lacks nominal-span predecessor evidence: "+read.id);
			}
			// A seed-terminal count cannot validate a candidate later in the added suffix.
			if(retained!=seedLength || read.length()==seedLength || actual>2 ||
					!rejectBackwardStep(context, actual, other, 2, 3, 1)){return false;}
			evidence[0]=context; evidence[1]=actual; evidence[2]=other;
			evidence[3]=otherBase; evidence[4]=tableK; evidence[5]=seedLength;
			return true;
		}

		/** Extra rows indicate a mismatched query population, not harmless unused evidence. */
		void requireEnd(){
			if(file.nextLine()!=null){throw new IllegalArgumentException("Hollow candidates contain more rows than saved queries.");}
		}
		boolean close(){return file.close();}
		private final ByteFile file;
		private final LineParser1 lp=new LineParser1('\t');
		int solidDepth=-1;
	}

	/** Rejects malformed paths before counting a base or handing a batch to the writer. */
	private void validatePath(final Read read){
		assert(pathSeedLength>=pathSpans[0]) : "The largest predecessor window must lie entirely within the retained seed.";
		if(read.mate!=null || read.length()<pathSeedLength || read.id==null || read.id.indexOf('\t')>=0){
			throw new IllegalArgumentException("Path veto requires an unpaired record at least pathseedlen long and a tab-free ID: "+read.id);
		}
		final int contextSpan=pathBackStop || pathDropStop || pathRiseStop || pathLowBranchStop || pathSpacingAbstain ? Math.max(k+1, pathSpans[0]) : pathSpans[0];
		for(int i=pathSeedLength-contextSpan; i<read.length(); i++){
			if(AminoAcid.baseToNumber[read.bases[i]]<0){throw new IllegalArgumentException("Ambiguous path query: "+read.id);}
		}
	}

	/** Confirms a weak Hollow backward seed-start candidate using its complete Solid word.
	 * This is a separate experimental conjunction, not a general coverage requirement.
	 * Strong Solid support retains even a Hollow candidate; no seed base is ever edited. */
	private int solidBackwardLimit(final Read read, final Kmer word, final int[] evidence, final CandidateCursor candidates,
			final Kmer shortWord, final Kmer shortNeighbor, final int[] fork){
		assert(pathSolidBackStop && !tables.hollow && candidates!=null && pathSeedLength>k) :
				"The exclusive seed-start diagnostic requires a Solid table and an ordered Hollow candidate stream.";
		Arrays.fill(evidence, -1);
		if(pathBackFork){Arrays.fill(fork, -1);}
		evidence[4]=k;
		if(!candidates.consume(read, pathSeedLength, k, evidence)){return read.length();}
		word.clear();
		for(int i=pathSeedLength-k; i<pathSeedLength; i++){word.addRight(read.bases[i]);}
		candidates.solidDepth=depth(word);
		if(!rejectSolidBackwardSeed(evidence[0], evidence[1], evidence[2], candidates.solidDepth)){return read.length();}
		if(pathBackFork){
			final int retained=seedForkLimit(read, word, shortWord, shortNeighbor, fork);
			if(retained<read.length()){evidence[6]=SOLID_SEED_SHORT_FORK;}
			return retained;
		}
		evidence[6]=SOLID_WEAK_BACKWARD;
		return pathSeedLength;
	}

	/** Confirms a weak seed-start candidate with a supported short-key fork.
	 * The full word is already positioned at the seed boundary. Its recovery ends
	 * this extra check permanently; absence alone never truncates. The k-64 bound
	 * follows the ignored interior width and cannot cross the saved read's end. */
	private int seedForkLimit(final Read read, final Kmer word, final Kmer shortWord,
			final Kmer neighbor, final int[] fork){
		assert(pathBackFork && !word.hollow && word.kbig==k && word.len>=k &&
				shortTables!=null && !shortTables.hollow && shortTables.kbig==64 && fork.length==5) :
				"Seed candidate validation supplies a complete SolidK context and separate Solid64 evidence; no Hollow counts may enter this fork test.";
		shortWord.clear();
		for(int i=pathSeedLength-64; i<pathSeedLength; i++){shortWord.addRight(read.bases[i]);}
		final int limit=seedForkHorizon(pathSeedLength, read.length(), k);
		for(int end=pathSeedLength; end<limit; end++){
			if(depth(word)>=2){return read.length();}
			final int chosenBase=AminoAcid.baseToNumber[read.bases[end]];
			final int context=Math.max(0, shortTables.getCount(shortWord));
			if(context>=3){
				neighbor.setFrom(shortWord); neighbor.addRightNumeric(chosenBase);
				final int chosen=Math.max(0, shortTables.getCount(neighbor));
				if(chosen>=3){
					int other=0, otherBase=-1;
					for(int base=0; base<4; base++){
						if(base==chosenBase){continue;}
						neighbor.setFrom(shortWord); neighbor.addRightNumeric(base);
						final int count=Math.max(0, shortTables.getCount(neighbor));
						if(count>other){other=count; otherBase=base;}
					}
					if(rejectSeedFork(context, chosen, other)){
						fork[0]=context; fork[1]=chosen; fork[2]=other;
						fork[3]=otherBase; fork[4]=end;
						return end;
					}
				}
			}
			word.addRightNumeric(chosenBase);
			shortWord.addRightNumeric(chosenBase);
		}
		return read.length();
	}

	/** Bounds uncertainty to the Hollow interior width without overflowing long reads. */
	public static int seedForkHorizon(final int seedLength, final int readLength, final int tableK){
		assert(tableK>64 && seedLength>tableK && readLength>=seedLength) :
				"The candidate's full seed context must fit, and the uncertainty interval may inspect added bases only.";
		return (int)Math.min((long)readLength, (long)seedLength+tableK-64);
	}

	/** A short fork matters only inside an independently established weak-seed interval.
	 * Both branches need repeated observations; equality and stronger alternatives
	 * stop, but a missing short context or chosen key must abstain. */
	public static boolean rejectSeedFork(final int context, final int chosen, final int other){
		assert(context>=0 && chosen>=0 && other>=0) :
				"Solid64 counts normalize absent sentinels before the overflow-safe supported-fork comparison.";
		return context>=3 && chosen>=3 && other>=3 && 2L*other>=chosen;
	}

	/** Stops where a supported nominal context contradicts an unsupported saved next base.
	 * Missing coverage abstains; the next context can regain support later in the path.
	 * Solid and Hollow modes use their own support floors and identical rolling logic.
	 * Optional nominal Hollow rises use the existing rise predicate before a positive
	 * saved choice. Optional source qualification requires that choice to remain weak
	 * at its original span, rather than rejecting legitimate seed-error recovery.
	 * Optional weak-choice dominance admits count1 choices only, retaining the same
	 * context/alternative support floors and any requested source qualification.
	 * Optional Solid ambiguity stops require both choices to meet the Solid floor;
	 * missing coverage still abstains, and no replacement base is selected.
	 * Optional spacing abstention protects a supported saved choice at k-1 or k+1
	 * observed spacing, queried against this same immutable nominal table.
	 * Optional read-only short support additionally requires a weak contained short word
	 * or weak proposed short successor. Declined cuts still advance every rolling state.
	 * Only original input evidence is counted, and the seed prefix is never changed. */
	private int contradictionLimit(final Read read, final Kmer word, final Kmer neighbor, final int[] evidence,
			final Kmer sourceWord, final int[] sourceEvidence, final Kmer spacingWord, final ShortKmerPathSupport shortSupport){
		assert(((pathSolidStop && !tables.hollow) || (pathNominalStop && tables.hollow)) &&
				querySpan==k && pathSeedLength>=k && read.length()>=pathSeedLength && evidence.length==7) :
				"Exclusive contradiction modes require matching table/query geometry, a complete seed window and seven evidence slots.";
		final int support=pathNominalStop ? pathSupport : pathSolidSupport;
		Arrays.fill(evidence, -1);
		if(sourceEvidence!=null){sourceEvidence[0]=-1;}
		evidence[4]=k;
		word.clear();
		for(int i=pathSeedLength-k; i<pathSeedLength; i++){word.addRight(read.bases[i]);}
		if(shortSupport!=null){shortSupport.reset(read.bases, pathSeedLength);}
		int context=depth(word);
		for(int end=pathSeedLength; end<read.length(); end++){
			final int chosenBase=AminoAcid.baseToNumber[read.bases[end]];
			neighbor.setFrom(word); neighbor.addRightNumeric(chosenBase);
			final int chosen=depth(neighbor);
			if(pathNominalRise && risePathStep(context, chosen, pathRiseRatio, pathSupport, pathRiseMin)){
				final int sourceChosen=pathRiseSource ? sourceChoiceCount(read, end, sourceWord) : -1;
				if(!pathRiseSource || weakSourceChoice(sourceChosen, pathSupport)){
					evidence[0]=context; evidence[1]=chosen; evidence[5]=end;
					evidence[6]=NOMINAL_HOLLOW_RISE;
					if(pathRiseSource){sourceEvidence[0]=sourceChosen;}
					return end;
				}
			}
			if(context>=support && (chosen==0 || (pathNominalWeak && chosen==1) ||
					(pathSolidBranchStop && chosen>=support))){
				int otherCount=0, otherBase=-1;
				for(int base=0; base<4; base++){
					if(base==chosenBase){continue;}
					neighbor.setFrom(word); neighbor.addRightNumeric(base);
					final int count=depth(neighbor);
					if(count>otherCount){otherCount=count; otherBase=base;}
				}
				final boolean rejected=pathSolidStop ? rejectSolidStep(context, chosen, otherCount, support, pathRatio, pathSolidBranchStop) :
						rejectUnsupportedStep(context, chosen, otherCount, support) ||
						(pathNominalWeak && rejectNominalWeakStep(context, chosen, otherCount, pathRatio, support));
				if(rejected){
					final int sourceChosen=sourceTables==null ? -1 : sourceChoiceCount(read, end, sourceWord);
					if((sourceTables==null || sourceChoiceOutvoted(otherCount, sourceChosen)) &&
							(!pathSpacingAbstain || !supportedSpacingChoice(read, end, spacingWord)) &&
							(shortSupport==null || shortSupport.allowsCut(chosenBase))){
						evidence[0]=context; evidence[1]=chosen; evidence[2]=otherCount;
						evidence[3]=otherBase; evidence[5]=end;
						evidence[6]=pathSolidStop && chosen>0 ? SOLID_AMBIGUITY : chosen==1 ? NOMINAL_HOLLOW_WEAK :
							(pathNominalStop ? NOMINAL_HOLLOW_CONTRADICTION : SOLID_CONTRADICTION);
						if(sourceEvidence!=null){sourceEvidence[0]=sourceChosen;}
						return end;
					}
				}
				// Restore the actual next word after using neighbor scratch for alternatives.
				neighbor.setFrom(word); neighbor.addRightNumeric(chosenBase);
			}
			word.setFrom(neighbor);
			context=chosen;
			if(shortSupport!=null){shortSupport.advance(chosenBase);}
		}
		return read.length();
	}

	/** Tests indel spacing only after a nominal contradiction would otherwise cut.
	 * Each observed context and its saved next key must independently meet the same
	 * support floor. Missing evidence offers no protection. Scratch is rebuilt for
	 * each span; the caller's nominal rolling word and neighbor remain untouched.
	 * No seed base, replacement choice, extra table or per-position allocation is needed. */
	private boolean supportedSpacingChoice(final Read read, final int end, final Kmer word){
		assert(pathSpacingAbstain && tables.hollow && word!=null && word.hollow && word.kbig==k &&
				k>=67 && end>=pathSeedLength && pathSeedLength>k && end<read.length()) :
				"Constructor and path validation ensure both alternate contexts fit before the saved next base and use the nominal table's Hollow geometry.";
		for(int span=k-1; span<=k+1; span+=2){
			setHollowQuery(word, read.bases, end-span, span);
			if(depth(word)<pathSupport){continue;}
			word.addRight(read.bases[end]);
			if(depth(word)>=pathSupport){return true;}
		}
		return false;
	}

	/** Tests count1 nominal choices with the existing alternative-dominance ratio.
	 * Zero-choice contradictions keep their old predicate; count2+ choices abstain.
	 * This experimental veto truncates a suffix, never substitutes the alternative. */
	public static boolean rejectNominalWeakStep(final int context, final int chosen, final int other,
			final int ratio, final int support){
		assert(context>=0 && chosen>=0 && other>=0 && ratio>=2 && support>=2) :
				"Nominal weak-choice counts are normalized and need independent repeated context/alternative support before a dominance comparison.";
		return chosen==1 && rejectPathStep(context, chosen, other, ratio, support);
	}

	/** Counts a saved choice at its original shorter Hollow span only for a candidate veto.
	 * Rare candidate queries rebuild one scratch word; ordinary positions incur no
	 * second-table lookup. A declined candidate does not end inspection of later bases. */
	private int sourceChoiceCount(final Read read, final int end, final Kmer word){
		assert(pathNominalStop && sourceTables!=null && word!=null && word.hollow &&
				word.kbig==pathSourceK && pathSourceK<k && end>=pathSeedLength && end<read.length()) :
				"Source support must use the shorter original-input geometry at the same saved next-base position.";
		word.clear();
		for(int i=end-pathSourceK+1; i<=end; i++){word.addRight(read.bases[i]);}
		return Math.max(0, sourceTables.getCount(word));
	}

	/** A longer-span contradiction must strictly outweigh the supported source choice.
	 * Equal counts can represent different repeat copies, as in the mixed-error g0490
	 * false stop; a tie does not justify discarding the original-span path.
	 * This raw-count comparison remains experimental, not a calibrated probability. */
	public static boolean sourceChoiceOutvoted(final int alternative, final int sourceChosen){
		assert(alternative>=0 && sourceChosen>=0) :
				"Both immutable tables normalize absent keys before comparing independent observations.";
		return sourceChosen>0 && alternative>sourceChosen;
	}

	/** A source-qualified rise requires positive but insufficient chosen-key support.
	 * Missing evidence abstains; reaching the existing support floor protects recovery.
	 * The source context is deliberately irrelevant to this next-base qualification. */
	public static boolean weakSourceChoice(final int chosen, final int support){
		assert(chosen>=0 && support>=2) :
				"Normalized original-span counts and a repeated-observation floor distinguish weak chosen evidence from absence or recovery.";
		return chosen>0 && chosen<support;
	}

	/** Stops at the first hard conflict; k+1 wins same-position ties.
	 * Optional rollback returns to the earliest ambiguous added-base choice in that
	 * same spacing, but ambiguity alone never truncates in rollback mode.
	 * The separate branch-stop experiment stops immediately at a supported ambiguity.
	 * Checkpoints cannot cross
	 * reads, spacings or the original seed boundary. Without rollback this remains
	 * the minimum retained length of independent single-spacing passes. */
	private int pathLimit(final Read read, final Kmer word, final Kmer neighbor, final int[] evidence, final int[] ambiguity){
		assert(read.length()>=pathSeedLength && evidence.length==7 && ambiguity.length==pathSpans.length) :
				"Validated reads and separate reusable evidence/checkpoints protect the seed and spacing-specific rollback history.";
		Arrays.fill(evidence, -1);
		Arrays.fill(ambiguity, -1);
		// Preserve the old single-span report even for untrimmed reads.
		if(!pathBothSpans){evidence[4]=querySpan;}
		for(int end=pathSeedLength; end<read.length(); end++){
			for(int index=0; index<pathSpans.length; index++){
				final int span=pathSpans[index];
				setHollowQuery(word, read.bases, end-span, span);
				final int context=depth(word);
				if(context<pathSupport){continue;}
				final int chosen=AminoAcid.baseToNumber[read.bases[end]];
				int chosenCount=0, otherCount=0, otherBase=-1;
				for(int base=0; base<4; base++){
					neighbor.setFrom(word); neighbor.addRightNumeric(base);
					final int count=depth(neighbor);
					if(base==chosen){chosenCount=count;}
					else if(count>otherCount){otherCount=count; otherBase=base;}
				}
				if(rejectPathStep(context, chosenCount, otherCount, pathRatio, pathSupport) ||
						(pathBranchStop && ambiguousPathStep(context, chosenCount, otherCount, pathRatio, pathSupport))){
					evidence[0]=context; evidence[1]=chosenCount; evidence[2]=otherCount;
					evidence[3]=otherBase; evidence[4]=span; evidence[5]=end;
					evidence[6]=ALTERNATE_SPAN;
					final int checkpoint=ambiguity[index];
					assert(checkpoint<0 || (checkpoint>=pathSeedLength && checkpoint<end)) :
							"A remembered ambiguity precedes this contradiction and must not discard original seed bases.";
					return checkpoint<0 ? end : checkpoint;
				}
				if(pathRollback && ambiguity[index]<0 && ambiguousPathStep(context, chosenCount, otherCount, pathRatio, pathSupport)){
					ambiguity[index]=end;
				}
			}
		}
		return read.length();
	}

	/** Checks nominal spacing before the existing retained boundary, without changing counts.
	 * A backward stop compares the observed predecessor against its alternatives, not
	 * the next base. Drop/rise stops compare the current window with the saved next window.
	 * Existing stops win equal-position ties; backward dominance precedes drop, rise,
	 * then the optional low-count forward-branch check. */
	private int nominalLimit(final Read read, final int retained, final Kmer word, final Kmer neighbor, final int[] evidence){
		assert(pathSeedLength>k && retained>=pathSeedLength && retained<=read.length()) :
				"Nominal checks require the preceding base and may only shorten an already validated extension suffix.";
		word.clear();
		for(int i=pathSeedLength-k; i<pathSeedLength; i++){word.addRight(read.bases[i]);}
		for(int end=pathSeedLength; end<retained; end++){
			final int context=depth(word);
			if(pathBackStop && context>=pathSupport){
				final int previous=AminoAcid.baseToNumber[read.bases[end-k-1]];
				int actualCount=0, otherCount=0, otherBase=-1;
				for(int base=0; base<4; base++){
					neighbor.setFrom(word); neighbor.addLeftNumeric(base);
					final int count=depth(neighbor);
					if(base==previous){actualCount=count;}
					else if(count>otherCount){otherCount=count; otherBase=base;}
				}
				if(rejectBackwardStep(context, actualCount, otherCount, pathBackRatio, pathSupport, pathBackActualSupport)){
					evidence[0]=context; evidence[1]=actualCount; evidence[2]=otherCount;
					evidence[3]=otherBase; evidence[4]=k; evidence[5]=end; evidence[6]=BACKWARD_DOMINANCE;
					return end;
				}
			}
			neighbor.setFrom(word); neighbor.addRight(read.bases[end]);
			if(pathDropStop || pathRiseStop || pathLowBranchStop){
				final int chosen=depth(neighbor);
				final boolean drop=pathDropStop && dropPathStep(context, chosen, pathDropRatio, pathDropContext, pathDropMax);
				final boolean rise=pathRiseStop && risePathStep(context, chosen, pathRiseRatio, pathSupport, pathRiseMin);
				if(drop || rise){
					evidence[0]=context; evidence[1]=chosen; evidence[2]=-1;
					evidence[3]=-1; evidence[4]=k; evidence[5]=end; evidence[6]=drop ? DEPTH_DROP : DEPTH_RISE;
					return end;
				}
				if(pathLowBranchStop && context>=pathSupport && chosen>=pathSupport){
					final int chosenBase=AminoAcid.baseToNumber[read.bases[end]];
					int otherCount=0, otherBase=-1;
					for(int base=0; base<4; base++){
						if(base==chosenBase){continue;}
						neighbor.setFrom(word); neighbor.addRightNumeric(base);
						final int count=depth(neighbor);
						if(count>otherCount){otherCount=count; otherBase=base;}
					}
					if(rejectLowBranchStep(context, chosen, otherCount, pathLowBranchRatio, pathSupport)){
						evidence[0]=context; evidence[1]=chosen; evidence[2]=otherCount;
						evidence[3]=otherBase; evidence[4]=k; evidence[5]=end; evidence[6]=LOW_FORWARD_BRANCH;
						return end;
					}
					// Neighbor scratch holds the last competitor; restore the saved path before advancing.
					neighbor.setFrom(word); neighbor.addRightNumeric(chosenBase);
				}
			}
			word.setFrom(neighbor);
		}
		return retained;
	}

	/** Requires two supported incoming paths before declaring backward dominance.
	 * A weak actual predecessor can be an existing seed error just outside the current
	 * window; that does not establish a hazard in the otherwise supported forward suffix. */
	public static boolean rejectBackwardStep(final int context, final int actual, final int other, final int ratio, final int floor){
		return rejectBackwardStep(context, actual, other, ratio, floor, floor);
	}

	/** Tests a separate actual-predecessor floor without weakening context or competitor support.
	 * Lowering this floor is experimental: it can reject both misjoined paths and valid
	 * extensions of seeds containing an error just outside the current window. */
	public static boolean rejectBackwardStep(final int context, final int actual, final int other,
			final int ratio, final int floor, final int actualFloor){
		assert(context>=0 && actual>=0 && other>=0 && ratio>=2 && floor>=1 && actualFloor>=1) :
				"Exact counts and positive independent floors protect the supported, overflow-safe backward comparison; absent actual keys must abstain.";
		return context>=floor && actual>=actualFloor && other>=floor && other>=(long)ratio*actual;
	}

	/** Stops only a sharp supported-to-low transition, not an already shallow path. */
	public static boolean dropPathStep(final int context, final int chosen, final int ratio, final int floor, final int maximum){
		assert(context>=0 && chosen>=0 && ratio>=2 && floor>=1 && maximum>=0) :
				"Exact counts and validated thresholds protect the overflow-safe depth-drop comparison.";
		return context>=floor && chosen<=maximum && context>=(long)ratio*chosen;
	}

	/** Detects abrupt entry into a much deeper endpoint-key class.
	 * A supported predecessor is required: weak seed-error keys can regain normal
	 * coverage when their error leaves an anchor. This is a conservative experimental
	 * stop, not proof of a repeat or permission to choose a different branch. */
	public static boolean risePathStep(final int context, final int chosen, final int ratio, final int support, final int minimum){
		assert(context>=0 && chosen>=0 && ratio>=2 && support>=1 && minimum>=1) :
				"Exact counts and validated positive thresholds protect the supported-predecessor, overflow-safe depth-rise comparison.";
		return context>=support && chosen>=minimum && chosen>=(long)ratio*context;
	}

	/** Tests only count2/3 nominal competitors without tightening tolerance of count1 noise.
	 * This narrow experiment is not a general junction classifier: deeper competitors
	 * remain the extender's responsibility. Ratio equality retains the saved choice. */
	public static boolean rejectLowBranchStep(final int context, final int chosen, final int other, final int ratio, final int floor){
		assert(context>=0 && chosen>=0 && other>=0 && ratio>=2 && floor>=1) :
				"Exact nonnegative counts and validated support protect the overflow-safe minor-branch comparison.";
		return context>=floor && chosen>=floor && other>=2 && other<=3 && chosen<(long)ratio*other;
	}

	/** Solid absence alone is not a veto: both the context and an alternative need evidence. */
	public static boolean rejectSolidStep(final int context, final int chosen, final int other, final int support){
		return rejectUnsupportedStep(context, chosen, other, support);
	}

	/** Adds symmetric supported ambiguity without imposing blanket Solid coverage.
	 * Off retains the old zero-choice contradiction. Ratio equality or a dominant
	 * positive choice on either side is not an ambiguity; neither selects a base. */
	public static boolean rejectSolidStep(final int context, final int chosen, final int other,
			final int support, final int ratio, final boolean branches){
		assert(support>=2 && ratio>=1 && (!branches || ratio>=2)) :
				"Solid vetoes require repeated observations; enabled ambiguity requires a meaningful strict dominance ratio.";
		return rejectSolidStep(context, chosen, other, support) ||
				(branches && ambiguousPathStep(context, chosen, other, ratio, support));
	}

	/** Requires repeated context and alternative evidence; any observed saved choice abstains.
	 * Counts must all come from the same immutable table and nominal query geometry. */
	public static boolean rejectUnsupportedStep(final int context, final int chosen, final int other, final int support){
		assert(context>=0 && chosen>=0 && other>=0 && support>=2) :
				"Normalized counts and a multi-observation floor prevent absence or one observation from establishing a contradiction.";
		return context>=support && chosen==0 && other>=support;
	}

	/** Combines fixed weak-predecessor Hollow evidence with insufficient full-Solid seed support. */
	public static boolean rejectSolidBackwardSeed(final int context, final int actual, final int other, final int solid){
		assert(solid>=0) : "The exact Solid table's absent sentinel is normalized before testing seed support.";
		return rejectBackwardStep(context, actual, other, 2, 3, 1) && actual<=2 && solid<2;
	}

	/** A supported alternative spacing may veto, but never replace, the saved next base. */
	public static boolean rejectPathStep(final int context, final int chosen, final int other, final int ratio, final int floor){
		assert(context>=0 && chosen>=0 && other>=0 && ratio>=1 && floor>=1) :
				"Exact counts and validated thresholds must remain nonnegative before the overflow-safe ratio comparison.";
		return context>=floor && other>=floor && other>=(long)ratio*chosen;
	}

	/** Remembers supported competing choices only; neither side may dominate by the veto ratio. */
	public static boolean ambiguousPathStep(final int context, final int chosen, final int other, final int ratio, final int floor){
		assert(context>=0 && chosen>=0 && other>=0 && ratio>=1 && floor>=1) :
				"Exact nonnegative counts and validated thresholds protect the symmetric, overflow-safe ambiguity comparison.";
		return context>=floor && chosen>=floor && other>=floor &&
				chosen<(long)ratio*other && other<(long)ratio*chosen;
	}

	/** Reuses observed endpoint anchors at another table span, without changing their bases.
	 * The extra base on each arm preserves one-step neighbor queries. All synthetic A bases
	 * remain strictly inside the ignored interior, even after either one-base shift. */
	public static void setHollowQuery(final Kmer word, final byte[] bases, final int start, final int span){
		if(!word.hollow || word.kbig<66 || span<66 || start<0 || start>bases.length-span){
			throw new IllegalArgumentException("Alternate-span queries require valid bounds and two nonoverlapping33bp contexts.");
		}
		assert(word.hollow && word.kbig>=66) : "Synthetic interior bases must not enter the central or one-step neighbor fingerprints.";
		word.clear();
		for(int pos=start; pos<start+33; pos++){word.addRight(bases[pos]);}
		for(int pos=33; pos<word.kbig-33; pos++){word.addRightNumeric(0);}
		for(int pos=start+span-33; pos<start+span; pos++){word.addRight(bases[pos]);}
		assert(word.len==word.kbig) : "Every alternate-span query must fully populate its table-sized Kmer.";
	}

	/** Normalizes the table's absent-key sentinel to a zero observation count. */
	private int depth(final Kmer word){
		assert(word.len>=k) : "Exact depth queries cannot use incomplete reference words.";
		return Math.max(0, tables.getCount(word));
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String reference, scaffold, output, queries, sequenceOutput, selfSeeds;
	/** Explicit original FASTQ for the optional read-only Solid protection table. */
	private String pathShortReads;
	/** Exact original-read support geometry; the default preserves the Solid64 diagnostic. */
	private int pathShortK=64;
	/** Two-pass singleton removal preserves the <=1 class; filtered zero is not exact absence. */
	private boolean pathShortPrefilter=false;
	private final ArrayList<int[]> regions=new ArrayList<int[]>();
	private final KmerTableSetU tables;
	private final KmerTableSetU shortTables;
	private final KmerTableSetU sourceTables;
	private final int k;
	private int querySpan=0;
	private int queryStart=0;
	private int pathSeedLength=0, pathRatio=4, pathSupport=3;
	private boolean pathBothSpans=false, pathRollback=false, pathBranchStop=false;
	private boolean pathSelfStop=false;
	private boolean pathBackStop=false, pathDropStop=false, pathRiseStop=false;
	private boolean pathLowBranchStop=false;
	private int pathLowBranchRatio=32;
	private boolean pathSolidStop=false;
	/** Default-off supported ambiguity modifier for the read-only Solid postpass. */
	private boolean pathSolidBranchStop=false;
	private boolean pathNominalStop=false;
	/** Separately enables ratio-dominated count1 choices in the nominal Hollow scan. */
	private boolean pathNominalWeak=false;
	/** Protects supported alternate-spacing choices only in the zero-choice nominal diagnostic. */
	private boolean pathSpacingAbstain=false;
	/** Adds the existing supported-rise predicate to the exclusive nominal Hollow scan. */
	private boolean pathNominalRise=false;
	/** Qualifies nominal rises using the already-loaded original-span chosen key. */
	private boolean pathRiseSource=false;
	/** Zero disables the extra immutable Hollow table and preserves the original diagnostic. */
	private int pathSourceK=0;
	private boolean pathShortStop=false;
	private boolean pathShortDrop=false;
	private boolean pathSolidBackStop=false;
	private boolean pathBackFork=false;
	private String pathCandidates;
	private int pathSolidSupport=2;
	private int pathBackRatio=2, pathDropRatio=16, pathDropContext=32, pathDropMax=2;
	/** Zero inherits pathSupport after parsing, preserving the original coupled default. */
	private int pathBackActualSupport=0;
	private int pathRiseRatio=16, pathRiseMin=32;
	private final int[] pathSpans;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final int ALTERNATE_SPAN=0, SELF_ONLY=1, BACKWARD_DOMINANCE=2, DEPTH_DROP=3, DEPTH_RISE=4;
	private static final int LOW_FORWARD_BRANCH=5;
	private static final int SOLID_CONTRADICTION=6;
	private static final int SOLID_WEAK_BACKWARD=7;
	private static final int SOLID_SEED_SHORT_FORK=8;
	private static final int NOMINAL_HOLLOW_CONTRADICTION=9;
	private static final int NOMINAL_HOLLOW_RISE=10;
	private static final int NOMINAL_HOLLOW_WEAK=11;
	private static final int SOLID_AMBIGUITY=12;
	private static final String[] STOP_REASONS={"alternate_span", "self_only", "backward_dominance", "depth_drop", "depth_rise", "low_forward_branch", "solid_contradiction", "solid_weak_backward_seed", "solid_seed_short_fork", "nominal_hollow_contradiction", "nominal_hollow_rise", "nominal_hollow_weak_choice", "solid_ambiguity"};
	private static final String CANDIDATE_SCHEMA="id\tseed_length\toriginal_length\tretained_length\tk\tquery_span\tcontext_depth\tchosen_count\tother_count\tother_base\tself_count\texternal_count\tstop_reason";
}
