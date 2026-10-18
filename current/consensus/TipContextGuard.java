package consensus;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parser;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Bounded, experimental substitution guard for one selected contig end.
 * Original support comes from an unpadded alignment model; additional support
 * comes from a separately mapped context model. Neither model is changed.
 * Unique means one compatible assembly END overlap within the stated bounds,
 * not genomic uniqueness. A one-read true allele remains indistinguishable
 * from a one-read error. This class is not called by ordinary assembly/polishing.
 * @author Fischl
 */
public final class TipContextGuard {

	/** Uses native argument expansion and closes a redirected diagnostic stream. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, TipContextGuard.class, false);
		try{new TipContextGuard(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires separate immutable models and an assembly containing the target. */
	private TipContextGuard(final String[] args){
		assert(args!=null) : "PreParser supplies expanded arguments before path validation.";
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("originalmodel")){originalModel=b;}
			else if(a.equals("contextmodel")){contextModel=b;}
			else if(a.equals("originalkmers")){originalKmers=b;}
			else if(a.equals("assembly")){assembly=b;}
			else if(a.equals("out")){output=b;}
			else if(a.equals("report")){report=b;}
			else if(a.equals("offset")){offset=Integer.parseInt(b);}
			else if(a.equals("side")){side=b;}
			else if(a.equals("minoverlap")){minOverlap=Integer.parseInt(b);}
			else if(a.equals("maxoverlap")){maxOverlap=Integer.parseInt(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		check(originalModel!=null && contextModel!=null && assembly!=null && output!=null && report!=null,
				"Required: originalmodel= contextmodel= assembly= out= report= offset=");
		check(offset>=0 && minOverlap>=2*ANCHOR+1 && maxOverlap>=minOverlap,
				"Nonnegative context offset and overlap bounds >=33 are required.");
		check("left".equals(side) || "right".equals(side), "side must be left or right.");
		check(!parser.append, "Appending would mix independent guard experiments.");
		check(Tools.testInputFiles(false, true, originalModel, contextModel, assembly, originalKmers), "Unreadable guard inputs.");
		check(Tools.testOutputFiles(false, false, false, output, report), "Guard outputs must be new files.");
		check(Tools.testForDuplicateFiles(true, originalModel, contextModel, assembly, originalKmers, output, report),
				"Guard outputs must not alias their evidence inputs.");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Discovers candidates before making any change; ambiguous evidence causes abstention. */
	private void process(){
		final BaseGraph original=readModel(originalModel), context=readModel(contextModel);
		final KmerEvidence kmers=originalKmers==null ? null : readKmers(originalKmers, original);
		check(original.pad==0 && context.pad==0, "Only unpadded ConsensusMaker models have the stated coordinates.");
		check((long)offset+original.original.length<=context.original.length, "Target span exceeds context model.");
		for(int i=0; i<original.original.length; i++){
			check(original.original[i]==context.original[offset+i], "Context changed the original target at "+i);
		}
		final boolean left=side.equals("left");
		final byte[] target=left ? reverse(original.original) : original.original;
		final CandidateSet candidates=discover(target, original.name, left, original.original);
		final Candidate candidate=candidates.count==1 ? candidates.only : null;
		final boolean contextOK=candidate!=null && matchesContext(candidate, context.original, offset, target.length, left);
		final byte[] result=original.original.clone(), scratch=new byte[2];
		int proposed=0, accepted=0;
		final ByteStreamWriter audit=new ByteStreamWriter(report, false, false, false);
		audit.start();
		audit.print("#end_candidates="+candidates.count+"\tmin_overlap="+minOverlap+"\tmax_overlap="+maxOverlap+"\n");
		audit.print("#unique_within_bounds_only; original support excludes the reference pseudocount\n");
		audit.print("position1\toriginal\tproposal\toriginal_observations\tcontext_proposal_observations\t"
				+(kmers==null ? "" : "original_kmer_support\t")+"decision\n");
		final int oldDepth=ConsensusObject.minDepth;
		final float oldMaf=ConsensusObject.MAF_sub;
		final boolean oldOnlyNs=ConsensusObject.onlyConvertNs;
		try{
			ConsensusObject.minDepth=2;
			ConsensusObject.MAF_sub=0.5f;
			ConsensusObject.onlyConvertNs=false;
			final ByteBuilder row=new ByteBuilder();
			for(int pos=0; pos<result.length; pos++){
				final byte from=original.original[pos], to=context.ref[offset+pos].consensus(scratch);
				if(to==from){continue;}
				proposed++;
				final int originalCount=observations(original.ref[pos], from);
				final int contextCount=observations(context.ref[offset+pos], to);
				final int orientedPos=left ? result.length-1-pos : pos;
				final byte orientedTo=left ? AminoAcid.baseToComplementExtended[to] : to;
				String reason=decision(originalCount, contextCount, orientedPos, orientedTo,
						target.length, candidates, contextOK);
				final int kmerSupport=kmers==null ? -1 : kmers.support(pos);
				if(reason.equals("accept") && kmers!=null){
					if(kmerSupport>=2){reason="preserve_original_kmer_support";}
					else if(kmerSupport<0){reason="missing_original_kmer_evidence";}
				}
				if(reason.equals("accept")){result[pos]=to; accepted++;}
				row.clear();
				row.append(pos+1).tab().append((char)from).tab().append((char)to).tab()
						.append(originalCount).tab().append(contextCount).tab();
				if(kmers!=null){row.append(kmerSupport).tab();}
				row.append(reason).nl();
				audit.print(row);
			}
		}finally{
			ConsensusObject.minDepth=oldDepth;
			ConsensusObject.MAF_sub=oldMaf;
			ConsensusObject.onlyConvertNs=oldOnlyNs;
			check(!audit.poisonAndWait(), "Failed to finish guard decision report.");
		}
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
		writer.start();
		writer.print(new ByteBuilder().append('>').append(original.name).nl().append(result).nl());
		check(!writer.poisonAndWait(), "Failed to finish guarded FASTA.");
		System.err.println("TIP_CONTEXT_GUARD candidates="+candidates.count+" proposed="+proposed+" accepted="+accepted);
	}

	/** Streams whole contigs once; only a bounded candidate sequence is retained. */
	private CandidateSet discover(final byte[] target, final String targetName, final boolean left, final byte[] storedTarget){
		assert(target.length==storedTarget.length) : "Orientation must preserve target coordinates and length.";
		FastaReadInputStream.SPLIT_READS=false;
		final FileFormat ff=FileFormat.testInput(assembly, FileFormat.FASTA, null, true, true);
		check(ff.fasta(), "The candidate universe must be a FASTA assembly, not reads or alignments.");
		final Streamer streamer=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 1);
		final CandidateSet set=new CandidateSet();
		int targets=0;
		streamer.start();
		try{
			for(ListNum<Read> batch=streamer.nextList(); batch!=null; batch=streamer.nextList()){
				for(Read read : batch){
					check(read.mate==null && read.bases!=null, "Assembly records must be unpaired concrete contigs.");
					if(read.id.equals(targetName)){
						check(Arrays.equals(read.bases, storedTarget), "Assembly target differs from original model.");
						targets++;
						continue;
					}
					// Both donor orientations are scanned regardless of the stored target strand.
					scan(target, read.bases, minOverlap, maxOverlap, set);
					scan(target, reverse(read.bases), minOverlap, maxOverlap, set);
				}
			}
		}finally{streamer.close();}
		check(!streamer.errorState(), "Assembly streaming failed; candidate completeness is unknown.");
		check(targets==1, "Exactly one assembly record must match the original model name; found "+targets);
		return set;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts physical oriented placements, including identical competing copies. */
	static void scan(final byte[] target, final byte[] donor, final int min, final int max, final CandidateSet set){
		assert(target!=null && donor!=null && set!=null && min>=2*ANCHOR+1 && max>=min)
				: "Bounded end matching needs two anchors and initialized candidate accounting.";
		final int limit=Math.min(max, Math.min(target.length-1, donor.length-1));
		for(int overlap=min; overlap<=limit; overlap++){
			int mismatches=0, mismatch=-1;
			boolean defined=true;
			for(int i=0; i<overlap; i++){
				final byte a=target[target.length-overlap+i], b=donor[i];
				if(AminoAcid.baseToNumber[a]<0 || AminoAcid.baseToNumber[b]<0){defined=false; break;}
				if(a!=b){mismatch=i; if(++mismatches>1){break;}}
			}
			if(!defined || mismatches>1 || (mismatches==1 && (mismatch<ANCHOR || overlap-1-mismatch<ANCHOR))){continue;}
			set.count++;
			if(set.count==1){set.only=new Candidate(donor.clone(), overlap);}
		}
	}

	/** Requires every added context base on this side to come from the unique donor. */
	static boolean matchesContext(final Candidate candidate, final byte[] context, final int offset,
			final int length, final boolean left){
		assert(candidate!=null && offset>=0 && (long)offset+length<=context.length)
				: "Validated target coordinates are required before comparing added context.";
		final int extra=left ? offset : context.length-offset-length;
		if(extra<1 || extra>candidate.bases.length-candidate.overlap){return false;}
		for(int i=0; i<extra; i++){
			final byte base=left ? AminoAcid.baseToComplementExtended[context[offset-1-i]] : context[offset+length+i];
			if(base!=candidate.bases[candidate.overlap+i] || AminoAcid.baseToNumber[base]<0){return false;}
		}
		return true;
	}

	/** The guard cannot identify true one-read alleles; they intentionally remain a test counterexample. */
	static String decision(final int originalCount, final int contextCount, final int pos, final byte to,
			final int length, final CandidateSet candidates, final boolean contextOK){
		assert(originalCount>=0 && contextCount>=0 && pos>=0 && pos<length && candidates.count>=0)
				: "Evidence counts exclude pseudocounts; target coordinates must remain in bounds.";
		if(originalCount>=2){return "preserve_original_support";}
		if(candidates.count!=1){return "nonunique_or_missing_overlap";}
		if(!contextOK){return "context_not_from_unique_donor";}
		final Candidate candidate=candidates.only;
		final int donorPos=pos-(length-candidate.overlap);
		if(donorPos<0 || donorPos>=candidate.overlap){return "outside_overlap";}
		if(to!=candidate.bases[donorPos]){return "donor_disagrees";}
		if(contextCount<2){return "insufficient_alternate_support";}
		return "accept";
	}

	/** BaseGraph constructs one reference observation; it is not independent read support. */
	static int observations(final BaseNode node, final byte base){
		assert(node!=null && node.type==ConsensusObject.REF) : "Only reference nodes carry the known single pseudocount.";
		final int allele=AminoAcid.baseToNumber[base];
		if(allele<0){return 0;}
		final int count=node.acgtCount[allele]-(node.refBase==base ? 1 : 0);
		check(count>=0, "Missing reference pseudocount in original alignment model at "+node.rpos);
		return count;
	}

	/** Reads only trusted local ConsensusMaker output, never untrusted serialized objects. */
	private static BaseGraph readModel(final String path){
		final Object object=ReadWrite.readObject(path, false);
		check(object instanceof ArrayList<?>, "Expected ConsensusMaker graph list: "+path);
		final ArrayList<?> list=(ArrayList<?>)object;
		check(list.size()==1 && list.get(0) instanceof BaseGraph, "Guard requires one selected target per model: "+path);
		return (BaseGraph)list.get(0);
	}

	/** Validates exact query identity before using native read-only KmerWindowDepth observations. */
	private static KmerEvidence readKmers(final String path, final BaseGraph original){
		assert(original.pad==0) : "Native kmer-end coordinates refer to the unpadded original contig.";
		final KmerEvidence evidence=new KmerEvidence(original.original.length);
		final LineParser1 parser=new LineParser1('\t');
		final ByteBuilder word=new ByteBuilder(32);
		int rows=0;
		for(byte[] line : ByteFile.toLines(path)){
			if(line.length==0){continue;}
			parser.set(line);
			if(parser.termEquals("region", 0)){continue;}
			check(parser.terms()==17 && parser.parseInt(2)==32, "Expected native exact32mer evidence schema.");
			check(parser.termEquals(original.name, 1), "Kmer evidence names a different original contig.");
			final int start=parser.parseInt(3), end=parser.parseInt(4), depth=parser.parseInt(6);
			check(end>=32 && end<=original.original.length && start==end-31 && depth>=0,
					"Invalid original kmer coordinates/depth: "+start+"-"+end+" depth="+depth);
			word.clear();
			parser.appendTerm(word, 5);
			check(word.length()==32, "The evidence word must contain exactly32 bases.");
			for(int i=0; i<32; i++){
				check(word.array[i]==original.original[start-1+i], "Kmer query differs from original target at "+(start+i));
			}
			check(evidence.depths[end-1]<0, "Duplicate evidence for kmer end "+end);
			evidence.depths[end-1]=depth;
			rows++;
		}
		check(rows>0, "An empty kmer file cannot support correction.");
		return evidence;
	}

	/** Returns a reverse complement without mutating the caller's sequence. */
	static byte[] reverse(final byte[] bases){
		assert(bases!=null) : "An oriented contig must have concrete sequence storage.";
		final byte[] result=new byte[bases.length];
		for(int i=0; i<bases.length; i++){result[i]=AminoAcid.baseToComplementExtended[bases[bases.length-1-i]];}
		return result;
	}

	/** Keeps evidence/path failures loud even if assertions are accidentally disabled. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** One oriented suffix-to-prefix placement, not an assertion of biological origin. */
	static final class Candidate {
		Candidate(final byte[] bases_, final int overlap_){
			assert(bases_!=null && overlap_>0 && overlap_<bases_.length) : "A donor must extend beyond its overlap.";
			bases=bases_;
			overlap=overlap_;
		}
		final byte[] bases;
		final int overlap;
	}

	/** Multiplicity is preserved even when two donors have exactly the same sequence. */
	static final class CandidateSet {
		int count;
		Candidate only;
	}

	/** Mapping-independent veto only: repeated observations preserve an allele, never force a correction. */
	static final class KmerEvidence {
		KmerEvidence(final int length){
			check(length>=32, "The experimental kmer veto requires at least32 target bases.");
			depths=new int[length];
			Arrays.fill(depths, -1);
		}

		/** Returns strongest support, or -1 if a weak site's overlapping query set is incomplete. */
		int support(final int pos){
			assert(pos>=0 && pos<depths.length) : "A proposed substitution must retain original coordinates.";
			int maximum=0;
			boolean missing=false;
			final int last=pos+Math.min(31, depths.length-1-pos);
			for(int end=Math.max(31, pos); end<=last; end++){
				final int depth=depths[end];
				missing|=depth<0;
				maximum=Math.max(maximum, depth);
			}
			return maximum>=2 || !missing ? maximum : -1;
		}

		final int[] depths;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String originalModel, contextModel, originalKmers, assembly, output, report;
	private String side="right";
	private int offset=-1;
	private int minOverlap=60, maxOverlap=200;
	private static final int ANCHOR=16;
}
