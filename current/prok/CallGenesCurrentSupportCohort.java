package prok;

import java.util.ArrayList;

import dna.GeneticCode;
import fileIO.ByteStreamWriter;
import parse.Parse;
import prot.ProteinTailAssigner;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/**
 * Extracts the current stop-supported CallGenes HBM tail cohort for calibration.
 * This is a standalone diagnostic for gate design, not biological truth.
 * @author Yelan
 */
public final class CallGenesCurrentSupportCohort {

	public static void main(String[] args){
		ProteinSupport support=null;
		try{
			final Settings s=new Settings(args);
			support=ProteinSupport.load("core", s.reference, 0, 0, new String[]{s.out, s.outMap, s.outFamilies});
			new CallGenesCurrentSupportCohort(s, support.newWorker()).process();
		}catch(Throwable t){
			t.printStackTrace();
			System.exit(1);
		}finally{
			if(support!=null){support.close();}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private CallGenesCurrentSupportCohort(final Settings settings_, final ProteinSupport.Worker worker_){
		settings=settings_;
		worker=worker_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		final ArrayList<Read> reads=ReadInputStream.toReads(settings.in, fileIO.FileFormat.FA, -1);
		final ByteStreamWriter fasta=ByteStreamWriter.makeBSW(settings.out, settings.overwrite, false, true);
		final ByteStreamWriter mapping=ByteStreamWriter.makeBSW(settings.outMap, settings.overwrite, false, true);
		final ByteStreamWriter families=ByteStreamWriter.makeBSW(settings.outFamilies, settings.overwrite, false, true);
		mapping.println("query\tgenome\tcontig\tstart\tstop\tstrand\tfamily_id\trepresentative\ttail_length\tscore64\tidentity\tpaired_fraction\tcore_cov\tterminal_paired\tquery_suffix\treference_suffix\tquery_start\treference_start");
		families.println("query\tfamily_id");
		final ByteBuilder bb=new ByteBuilder();
		int emitted=0;
		for(Read read : reads){
			final ArrayList<Orf>[] frames=GeneCaller.makeOrfs(token(read.id), read.bases, settings.minLen, null);
			final byte[][] oriented={read.bases, AminoCopy.reverseComplement(read.bases)};
			for(int strand=0; strand<2; strand++){
				for(int frame=0; frame<3; frame++){
					final ArrayList<Orf> list=frames[frame+3*strand];
					if(list==null || list.isEmpty()){continue;}
					final ProteinSupport.StopEvidence[] evidence=worker.prepareStops(list, oriented[strand], settings.code);
					assert(evidence.length==list.size()) : "Stop evidence must preserve longest-ORF indexing";
					for(int i=0; i<evidence.length; i++){
						final ProteinSupport.StopEvidence stop=evidence[i];
						if(!stop.supported()){continue;}
						final Orf orf=list.get(i);
						assert(orf.flipped()==0 && orf.strand==strand) : "makeOrfs must return CDS coordinates on their native strand";
						final byte[] peptide=ProteinSupport.translateTail(oriented[strand], orf.start, orf.stop, settings.code, stop.length);
						assert(peptide.length==stop.length) : "Replayed peptide length must match ProteinSupport";
						bb.clear().append(settings.genome).append("_stop_");
						appendPadded(bb, emitted+1, 6);
						final String query=bb.toString();
						final int start=genomicStart(orf, read.length())+1, end=genomicStop(orf, read.length())+1;
						final ProteinTailAssigner.Hit hit=stop.hit;
						bb.clear().append('>').append(query)
							.append(" genome=").append(settings.genome)
							.append(" family=").append(hit.familyId)
							.append(" contig=").append(token(read.id))
							.append(" coords=").append(start).append('-').append(end)
							.append(" strand=").append(strand==0 ? '+' : '-')
							.append(" table=").append(settings.code.id()).nl()
							.append(peptide).nl();
						fasta.print(bb);
						families.print(query).tab().println(hit.familyId);
						bb.clear().append(query).tab().append(settings.genome).tab().append(token(read.id))
							.tab().append(start).tab().append(end).tab().append(strand==0 ? '+' : '-')
							.tab().append(hit.familyId).tab();
						AnnotationGff.escaped(bb, hit.representative);
						bb.tab().append(stop.length).tab().append(hit.score64)
							.tab().append(hit.identity, 8).tab().append(hit.pairedFraction, 8).tab().append(hit.coreCoverage, 8)
							.tab().append(hit.terminalPaired).tab().append(hit.querySuffix).tab().append(hit.referenceSuffix)
							.tab().append(hit.queryStart).tab().append(hit.referenceStart).nl();
						mapping.print(bb);
						emitted++;
					}
				}
			}
		}
		if(fasta.poisonAndWait()){throw new RuntimeException("FASTA write failed: "+settings.out);}
		if(mapping.poisonAndWait()){throw new RuntimeException("Mapping write failed: "+settings.outMap);}
		if(families.poisonAndWait()){throw new RuntimeException("Family map write failed: "+settings.outFamilies);}
		if(settings.expectedSearched>=0 && worker.searched!=settings.expectedSearched){
			throw new IllegalStateException("Expected searched="+settings.expectedSearched+" but observed "+worker.searched);
		}
		if(settings.expectedMatched>=0 && emitted!=settings.expectedMatched){
			throw new IllegalStateException("Expected matched="+settings.expectedMatched+" but observed "+emitted);
		}
		if(settings.expectedUnsupported>=0 && worker.unsupported!=settings.expectedUnsupported){
			throw new IllegalStateException("Expected unsupported="+settings.expectedUnsupported+" but observed "+worker.unsupported);
		}
		System.err.println("CALLGENES_CURRENT_SUPPORT_COHORT_PASS genome="+settings.genome+
			" searched="+worker.searched+" matched="+emitted+" unsupported="+worker.unsupported);
	}

	private static int genomicStart(final Orf orf, final int length){
		return orf.strand==0 ? orf.start : length-1-orf.stop;
	}

	private static int genomicStop(final Orf orf, final int length){
		return orf.strand==0 ? orf.stop : length-1-orf.start;
	}

	private static String token(final String id){
		final int idx=Tools.indexOfWhitespace(id);
		return idx<0 ? id : id.substring(0, idx);
	}

	private static void appendPadded(final ByteBuilder bb, final int number, final int digits){
		assert(number>=0 && digits>0) : "Cohort IDs use nonnegative numbers and a positive digit width";
		int places=1;
		for(int limit=10; places<digits && number>=limit; places++, limit*=10){}
		for(; places<digits; places++){bb.append('0');}
		bb.append(number);
	}

	private static final class Settings{
		Settings(final String[] args){
			for(String arg : args){
				final String[] split=arg.split("=", 2);
				final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
				if(a.equals("in")){in=b;}
				else if(a.equals("reference")){reference=b;}
				else if(a.equals("genome")){genome=b;}
				else if(a.equals("out")){out=b;}
				else if(a.equals("outmap")){outMap=b;}
				else if(a.equals("outfamilies")){outFamilies=b;}
				else if(a.equals("minlen")){minLen=Parse.parseIntKMG(b);}
				else if(a.equals("expectedsearched")){expectedSearched=Parse.parseIntKMG(b);}
				else if(a.equals("expectedmatched")){expectedMatched=Parse.parseIntKMG(b);}
				else if(a.equals("expectedunsupported")){expectedUnsupported=Parse.parseIntKMG(b);}
				else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
				else{throw new IllegalArgumentException("Unknown parameter "+arg);}
			}
			if(in==null || reference==null || genome==null || out==null || outMap==null || outFamilies==null){
				throw new IllegalArgumentException("Required: in= reference= genome= out= outmap= outfamilies=");
			}
			if(minLen<3){throw new IllegalArgumentException("minlen must be at least 3");}
		}
		String in, reference, genome, out, outMap, outFamilies;
		int minLen=80, expectedSearched=-1, expectedMatched=-1, expectedUnsupported=-1;
		boolean overwrite;
		final GeneticCode code=GeneticCode.forTable(11);
	}

	private static final class AminoCopy{
		static byte[] reverseComplement(final byte[] bases){
			final byte[] copy=bases.clone();
			dna.AminoAcid.reverseComplementBasesInPlace(copy);
			return copy;
		}
	}

	private final Settings settings;
	private final ProteinSupport.Worker worker;
}
