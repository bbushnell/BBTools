package prok;

import java.io.PrintWriter;
import java.util.ArrayList;

import dna.AminoAcid;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;

/**
 * Boundary-instrumentation capture driver for the r58/lsu scavenger families: runs the
 * production NcrnaScavenger pipeline over one genome with NcrnaBoundaryInstrumentSink armed,
 * emitting one TSV row per ACCEPTED locus (trimSucceeded/nnInvoked included) plus the final
 * accepted calls as GFF. This is the "L1" measurement of the rrnafallback boundary
 * investigation: it answers, per accepted call, whether the alignment-extent snap actually
 * ran (trimSucceeded=true) or the emitted span is the raw alignWindow span (=false, the
 * escape hatch documented at NcrnaScavenger.trimToAlignmentExtent).
 *
 * <p>Faithfulness contract: families are loaded through CallGenes.loadNcrnaResources() (the
 * production loader, same override statics as the callgenes r58consensus=/r58models=/... flags,
 * zero duplicated constants) and scavengers are constructed with the exact argument list
 * GeneCaller uses (GeneCaller.makeRnaOrfs, the ncrnaScavengers construction block), including
 * the same post-construction hbmPass/collapseFrac assignments and nothing else. Per
 * GeneCaller's own per-strand loop, each family gets a FRESH calledPos list per strand and the
 * minus pass runs on the reverse complement with results flipped afterward — reproduced here
 * verbatim. Families other than r58/lsu are deliberately NOT run: cross-family claiming is
 * independent by design (fresh calledPos per family per strand, see GeneCaller's comment in
 * makeRnaOrfs), so their absence cannot change r58/lsu scavenge-time behavior or captures.
 * Contigs are passed as-is (no case normalization), matching production CallGenes input.
 *
 * <p>Single-threaded by construction (one main-loop scavenge call at a time), satisfying the
 * sink contract's t=1 requirement.
 *
 * <p>SCOPE OF callsOut (Ganyu review, 2026-09-09): the GFF holds the DIRECT scavenger-returned
 * found lists, BEFORE any whole-GeneCaller final path arbitration — a driver-local record for
 * joining captures to spans, NOT claimed equivalent to a full callgenes run's emitted GFF.
 * The production model annotation NcrnaScavenger sets on each Orf is preserved untouched
 * (an earlier draft overwrote it with the family name — Ganyu catch; family identity is
 * derivable from the model-name prefix and is a first-class column in the captures TSV).
 *
 * <p>JOIN CONTRACT (Qiqi review): captures are PRE-NN and calls are post-NN with no shared ID.
 * In the current bootstrap state the boundary nets are structurally OFF (null templates), so
 * nn_invoked is always false and post-trim spans ARE the emitted spans — joins work on
 * (contig, family, abs coords). A future NN-on run must add an explicit shared locus ID before
 * trusting joins. As a loud driver-local postcondition, per-family capture counts must equal
 * per-family accepted-call counts (1:1 by construction: the sink fires once per accepted locus
 * inside trimToAlignmentExtent, and scavenge adds each accepted locus to the returned list);
 * a mismatch exits nonzero.
 *
 * <p>Output columns (captures TSV): genome, family, contig, strand (0/1, GeneCaller
 * convention), model_index, model_name, orig_wstart, orig_wstop, posttrim_start,
 * posttrim_stop (all four in the PER-STRAND frame the sink documents), abs_start, abs_stop
 * (contig-absolute, 0-based inclusive; strand=1 transformed by the self-inverse
 * scaflen-1-stop/scaflen-1-start formula from the sink javadoc), trim_succeeded, nn_invoked.
 *
 * @author G11
 */
public class NcrnaBoundaryCaptureDriver {

	public static void main(String[] args) throws Exception {
		String in=null, label=null, out=null, callsOut=null;
		for(String arg : args){
			final int eq=arg.indexOf('=');
			if(eq<1){throw new IllegalArgumentException("Expected key=value: "+arg);}
			final String key=arg.substring(0, eq).toLowerCase(), value=arg.substring(eq+1);
			if(key.equals("in")){in=value;}
			else if(key.equals("label")){label=value;}
			else if(key.equals("out")){out=value;}
			else if(key.equals("callsout")){callsOut=value;}
			else if(key.equals("r58kmers")){CallGenes.R58_KMERS_OVERRIDE=value;}
			else if(key.equals("r58consensus")){CallGenes.R58_CONSENSUS_OVERRIDE=value;}
			else if(key.equals("r58models")){CallGenes.R58_MODELS_OVERRIDE=value;}
			else if(key.equals("lsukmers")){CallGenes.LSU_KMERS_OVERRIDE=value;}
			else if(key.equals("lsuconsensus")){CallGenes.LSU_CONSENSUS_OVERRIDE=value;}
			else if(key.equals("lsumodels")){CallGenes.LSU_MODELS_OVERRIDE=value;}
			else{throw new IllegalArgumentException("Unknown parameter: "+key);}
		}
		if(in==null || label==null || out==null || callsOut==null){
			System.err.println("Usage: in=<genome.fna> label=<name> out=<captures.tsv> callsout=<calls.gff>");
			System.err.println("  r58kmers= r58consensus= r58models= lsukmers= lsuconsensus= lsumodels= (all six required)");
			System.exit(1);
		}

		CallGenes.NCRNA_FAMILIES_ENABLED=true;
		CallGenes.R58LSU_ENABLED=true;
		CallGenes.loadNcrnaResources();

		final ArrayList<NcrnaScavenger> scavengers=new ArrayList<>();
		final ArrayList<String> famNames=new ArrayList<>();
		for(NcrnaFamily fam : GeneCaller.ncrnaFamilies){
			if(!fam.name.equals("r58") && !fam.name.equals("lsu")){continue;}
			//Verbatim GeneCaller construction (makeRnaOrfs ncrnaScavengers block) — any drift
			//here invalidates the faithfulness contract in the class javadoc.
			NcrnaScavenger scavenger=new NcrnaScavenger(fam.library, fam.models, fam.modelNames,
				fam.kmerSet, fam.kLong, fam.minLen, fam.windowPad,
				fam.indexK, fam.indexTopN, fam.adaptive,
				fam.adaptFloor, fam.adaptTopFrac, fam.adaptQFrac, fam.fixedMinHits,
				fam.scoreA, fam.scoreB, fam.idPass, fam.idBorderline,
				fam.boundary5NetTemplate, fam.boundary3NetTemplate,
				fam.boundaryStartTable, fam.boundaryStopTable,
				fam.boundaryStartInside, fam.boundaryStartOutside, fam.boundaryStopInside, fam.boundaryStopOutside,
				fam.boundaryMeanLen, fam.boundaryStartOffsets, fam.boundaryStopOffsets);
			scavenger.hbmPass=fam.hbmPass;
			scavenger.collapseFrac=fam.collapseFrac;
			scavengers.add(scavenger);
			famNames.add(fam.name);
		}
		if(scavengers.size()!=2){
			throw new RuntimeException("Expected exactly r58+lsu families loaded, got "+famNames
				+" — check the six override paths (production requireR58LsuOverrides should have "
				+"caught a missing one; a wrong-but-present path fails later at resource load).");
		}
		//Explicit arming precondition, NOT gated on -ea (Qiqi review: setInstrumentSink's
		//alignment check is an assert; a -da invocation could silently emit calls without
		//captures). The production loader also validates this, but the driver's correctness
		//depends on it directly, so it fails loud here regardless of JVM flags.
		for(NcrnaFamily fam : GeneCaller.ncrnaFamilies){
			if(!fam.name.equals("r58") && !fam.name.equals("lsu")){continue;}
			if(fam.modelNames==null || fam.library==null || fam.modelNames.length!=fam.library.length){
				throw new RuntimeException("Family "+fam.name+" modelNames/library misaligned ("
					+(fam.modelNames==null ? "null" : fam.modelNames.length)+" vs "
					+(fam.library==null ? "null" : fam.library.length)+") — captures would be "
					+"silently skipped for some accepted loci (see setInstrumentSink's javadoc).");
			}
		}

		final String genomeLabel=label;
		try(PrintWriter capPw=new PrintWriter(out); PrintWriter callsPw=new PrintWriter(callsOut)){
			capPw.println("#genome\tfamily\tcontig\tstrand\tmodel_index\tmodel_name\torig_wstart\torig_wstop"
				+"\tposttrim_start\tposttrim_stop\tabs_start\tabs_stop\ttrim_succeeded\tnn_invoked");
			callsPw.println("##gff-version 3");
			//Mutable per-contig state read by the sinks; safe because everything below is one
			//thread (sink contract's t=1). scafLen is required for the strand=1 absolute transform.
			final int[] scafLen=new int[1];
			//Per-family counters for the capture==calls postcondition (Qiqi review): the sink
			//fires once per accepted locus inside trimToAlignmentExtent and scavenge() adds each
			//accepted locus to its returned list, so these must end equal; a mismatch means a
			//silent capture skip or an unaccounted call path and the run is invalid.
			final long[] captureCounts=new long[scavengers.size()];
			final long[] callCounts=new long[scavengers.size()];
			for(int fi=0; fi<scavengers.size(); fi++){
				final int fidx=fi;
				final String famName=famNames.get(fi);
				final NcrnaFamily fam=GeneCaller.ncrnaFamilies.stream()
					.filter(f -> f.name.equals(famName)).findFirst().get();
				scavengers.get(fi).setInstrumentSink(new NcrnaBoundaryInstrumentSink(){
					@Override
					public void capture(String contigName, int strand, int modelIndex,
							int origWStart, int origWStop, int postTrimStart, int postTrimStop,
							byte[] windowCopy, int windowCopyOffset, boolean trimSucceeded, boolean nnInvoked){
						//Self-inverse per-strand<->absolute transform from the sink javadoc
						//(PFeature.flip's formula): identity on strand 0.
						final int absStart=(strand==0 ? postTrimStart : scafLen[0]-1-postTrimStop);
						final int absStop=(strand==0 ? postTrimStop : scafLen[0]-1-postTrimStart);
						captureCounts[fidx]++;
						capPw.println(genomeLabel+"\t"+famName+"\t"+contigName+"\t"+strand
							+"\t"+modelIndex+"\t"+fam.modelNames[modelIndex]
							+"\t"+origWStart+"\t"+origWStop+"\t"+postTrimStart+"\t"+postTrimStop
							+"\t"+absStart+"\t"+absStop+"\t"+trimSucceeded+"\t"+nnInvoked);
					}
				});
			}

			final FileFormat ffin=FileFormat.testInput(in, FileFormat.FA, null, true, true);
			final ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, ffin, null);
			long contigs=0, calls=0;
			boolean streamError=false;
			try{
				cris.start();
				for(ListNum<Read> ln=cris.nextList(); ln!=null && ln.size()>0; ln=cris.nextList()){
					for(Read read : ln){
						if(read.bases==null){continue;}
						contigs++;
						scafLen[0]=read.length();
						final byte[] bases=read.bases;
						//rcBases computed once per contig, reused for all families' strand-1
						//passes — matching GeneCaller's cost profile and code shape.
						final byte[] rcBases=AminoAcid.reverseComplementBases(bases);
						for(int strand=0; strand<2; strand++){
							final byte[] strandBases=(strand==0 ? bases : rcBases);
							for(int fi=0; fi<scavengers.size(); fi++){
								ArrayList<int[]> calledPos=new ArrayList<>();
								ArrayList<Orf> found=scavengers.get(fi).scavenge(read.id, strandBases, strand, calledPos);
								if(found==null){continue;}
								if(strand==1){for(Orf orf : found){orf.flip();}}
								//Model annotation set by NcrnaScavenger is preserved untouched
								//(Ganyu review: an earlier draft overwrote it with the family
								//name, destroying the exact-model record production keeps).
								for(Orf orf : found){callsPw.println(orf.toGff());}
								callCounts[fi]+=found.size();
								calls+=found.size();
							}
						}
					}
					cris.returnList(ln);
				}
			}finally{streamError=ReadWrite.closeStream(cris);}
			//Read-side error visibility (Ganyu review, rev3): closeStream's boolean carries the
			//input stream's error state; ignoring it could pass a truncated/failed read off as a
			//complete genome. Checked HERE, after the finally, so a real exception from the try
			//body is never masked by a throw from within finally.
			if(streamError){
				throw new RuntimeException("Input stream error state on "+in
					+" — reads may be incomplete; captures/calls unreliable, do not consume.");
			}
			for(int fi=0; fi<scavengers.size(); fi++){
				System.err.println("family="+famNames.get(fi)+" captures="+captureCounts[fi]
					+" calls="+callCounts[fi]);
				if(captureCounts[fi]!=callCounts[fi]){
					throw new RuntimeException("POSTCONDITION FAIL: family "+famNames.get(fi)
						+" captures="+captureCounts[fi]+" != calls="+callCounts[fi]
						+" — a capture was silently skipped or a call path went unaccounted; "
						+"the capture TSV cannot be trusted as the accepted-denominator record.");
				}
			}
			System.err.println("NcrnaBoundaryCaptureDriver: genome="+genomeLabel+" contigs="+contigs
				+" calls="+calls);
			//PrintWriter swallows I/O errors (Qiqi review) — checkError() is the only visibility;
			//a buffered write failure here means partial/corrupt outputs, so no DONE sentinel.
			if(capPw.checkError() || callsPw.checkError()){
				throw new RuntimeException("WRITE ERROR on captures or calls output — files are "
					+"unreliable; do not consume.");
			}
		}
		System.err.println("CAPTURE_DRIVER_DONE "+genomeLabel);
	}
}
