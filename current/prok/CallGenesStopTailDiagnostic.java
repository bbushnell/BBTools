package prok;

import java.util.ArrayList;
import java.util.HashMap;

import dna.GeneticCode;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parse;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/**
 * Extracts actual longest-stop tails for a fixed set of old CallGenes CDS rows.
 * This is a standalone diagnostic for stop-tail calibration, not a caller flag.
 * @author Keqing, Yelan
 */
public final class CallGenesStopTailDiagnostic {

	public static void main(String[] args){
		try{
			final Settings s=new Settings(args);
			new CallGenesStopTailDiagnostic(s).process();
		}catch(Throwable t){t.printStackTrace(); System.exit(1);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private CallGenesStopTailDiagnostic(final Settings settings_){
		settings=settings_;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		final ArrayList<Target> targets=loadTargets(settings.map);
		final HashMap<String, ArrayList<Target>> byContig=group(targets);
		final ArrayList<Read> reads=ReadInputStream.toReads(settings.in, fileIO.FileFormat.FA, -1);
		final ByteStreamWriter fasta=ByteStreamWriter.makeBSW(settings.out, settings.overwrite, false, true);
		final ByteStreamWriter mapping=ByteStreamWriter.makeBSW(settings.outMap, settings.overwrite, false, true);
		final ByteStreamWriter families=ByteStreamWriter.makeBSW(settings.outFamilies, settings.overwrite, false, true);
		final ByteBuilder bb=new ByteBuilder();
		int matched=0, searchable=0, unsupported=0, absent=0;
		mapping.println("query\tcontig\ttarget_start\ttarget_stop\tstrand\tfamily_id\tlongest_start\tlongest_stop\toriented_start\toriented_stop\ttail_length\treason");
		families.println("query\tfamily_id");
		for(Read read : reads){
			final ArrayList<Target> list=byContig.get(token(read.id));
			if(list==null){continue;}
			final ArrayList<Orf>[] frames=GeneCaller.makeOrfs(token(read.id), read.bases, settings.minLen, null);
			final byte[][] oriented={read.bases, AminoCopy.reverseComplement(read.bases)};
			for(Target target : list){
				final Orf orf=find(frames, target, read.length());
				if(orf==null){
					mapping.print(bb.clear().append(target.id).tab().append(target.contig).tab().append(target.start1).tab().append(target.stop1)
						.tab().append(target.strand==0 ? '+' : '-').tab().append(target.family)
						.tab().append(-1).tab().append(-1).tab().append(-1).tab().append(-1).tab().append(0)
						.tab().append("no_longest_stop").nl());
					System.err.println("No longest ORF matched "+target.id+" at "+target.contig+
							":"+target.start1+"-"+target.stop1+"\n"+nearby(frames, target, read.length()));
					absent++;
					matched++;
					continue;
				}
				assert(orf.flipped()==0) : "makeOrfs must return CDS coordinates on their native strand";
				final int first=orf.start, last=orf.stop;
				final byte[] peptide=orf.stopTruncated ? new byte[0] : ProteinSupport.translateTail(oriented[orf.strand], first, last, settings.code, settings.cap);
				final String reason=orf.stopTruncated ? "synthetic_stop" : ProteinSupport.translationReason(peptide);
				if(reason==null){
					bb.clear().append('>').append(target.id).append(" family=").append(target.family)
						.append(" contig=").append(target.contig).append(" target=").append(target.start1).append('-').append(target.stop1)
						.append(" longest=").append(genomicStart(orf, read.length())+1).append('-').append(genomicStop(orf, read.length())+1)
						.append(" strand=").append(orf.strand==0 ? '+' : '-')
						.append(" table=").append(settings.code.id()).append(" cap=").append(settings.cap).nl()
						.append(peptide).nl();
					fasta.print(bb);
					families.print(target.id).tab().println(target.family);
					searchable++;
				}else{unsupported++;}
				mapping.print(bb.clear().append(target.id).tab().append(target.contig).tab().append(target.start1).tab().append(target.stop1)
					.tab().append(target.strand==0 ? '+' : '-').tab().append(target.family)
					.tab().append(genomicStart(orf, read.length())+1).tab().append(genomicStop(orf, read.length())+1)
					.tab().append(first).tab().append(last).tab().append(peptide.length)
					.tab().append(reason==null ? "searchable" : reason).nl());
				matched++;
			}
		}
		if(fasta.poisonAndWait()){throw new RuntimeException("FASTA write failed: "+settings.out);}
		if(mapping.poisonAndWait()){throw new RuntimeException("Mapping write failed: "+settings.outMap);}
		if(families.poisonAndWait()){throw new RuntimeException("Family map write failed: "+settings.outFamilies);}
		if(matched!=targets.size()){throw new IllegalStateException("Matched "+matched+" of "+targets.size()+" requested rows");}
		System.err.println("CALLGENES_STOP_TAIL_DIAGNOSTIC_PASS targets="+targets.size()+" searchable="+searchable+" unsupported="+unsupported+" absent="+absent);
	}

	private static Orf find(final ArrayList<Orf>[] frames, final Target target, final int length){
		Orf found=null;
		for(ArrayList<Orf> list : frames){
			if(list==null){continue;}
			for(Orf orf : list){
				if(!orf.scafName.equals(target.contig) || orf.strand!=target.strand){continue;}
				final int targetStop=target.strand==0 ? target.stop1-1 : length-target.start1;
				if(orf.stop!=targetStop){continue;}
				if(found!=null){throw new IllegalStateException("Duplicate longest ORF for "+target.id);}
				found=orf;
			}
		}
		return found;
	}

	private static String nearby(final ArrayList<Orf>[] frames, final Target target, final int length){
		final int targetStop=target.strand==0 ? target.stop1-1 : length-target.start1;
		final ByteBuilder bb=new ByteBuilder("Nearest same-strand longest ORFs:");
		int emitted=0;
		for(ArrayList<Orf> list : frames){
			if(list==null){continue;}
			for(Orf orf : list){
				if(!orf.scafName.equals(target.contig) || orf.strand!=target.strand){continue;}
				final int candidateStop=orf.stop;
				final int distance=Math.abs(candidateStop-targetStop);
				if(distance<=600){
					bb.nl().append(genomicStart(orf, length)+1).append('-').append(genomicStop(orf, length)+1)
						.append(" strand=").append(orf.strand==0 ? '+' : '-')
						.append(" stop_distance=").append(distance)
						.append(" flipped=").append(orf.flipped())
						.append(" stop_truncated=").append(orf.stopTruncated);
					emitted++;
				}
			}
		}
		if(emitted==0){bb.append(" none within 600 bp of biological stop");}
		return bb.toString();
	}

	private static int genomicStart(final Orf orf, final int length){
		return orf.strand==0 ? orf.start : length-1-orf.stop;
	}

	private static int genomicStop(final Orf orf, final int length){
		return orf.strand==0 ? orf.stop : length-1-orf.start;
	}

	private static ArrayList<Target> loadTargets(final String path){
		final ArrayList<Target> out=new ArrayList<Target>();
		final ByteFile bf=ByteFile.makeByteFile(path, false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				lp.set(line);
				if(lp.termEquals("query", 0)){continue;}
				if(lp.terms()<6){throw new IllegalArgumentException("Expected six tab-delimited target columns");}
				out.add(new Target(lp.parseString(0), lp.parseString(1), lp.parseInt(2), lp.parseInt(3), lp.parseString(4), lp.parseInt(5)));
			}
		}finally{
			if(bf.close()){throw new RuntimeException("Could not close target map: "+path);}
		}
		if(out.isEmpty()){throw new IllegalArgumentException("No diagnostic targets loaded");}
		return out;
	}

	private static HashMap<String, ArrayList<Target>> group(final ArrayList<Target> targets){
		final HashMap<String, ArrayList<Target>> map=new HashMap<String, ArrayList<Target>>();
		for(Target target : targets){
			ArrayList<Target> list=map.get(target.contig);
			if(list==null){map.put(target.contig, list=new ArrayList<Target>());}
			list.add(target);
		}
		return map;
	}

	private static String token(final String id){
		final int idx=Tools.indexOfWhitespace(id);
		return idx<0 ? id : id.substring(0, idx);
	}

	private static final class Settings{
		Settings(final String[] args){
			for(String arg : args){
				final String[] split=arg.split("=", 2);
				final String a=split[0].toLowerCase(), b=split.length>1 ? split[1] : null;
				if(a.equals("in")){in=b;}
				else if(a.equals("map")){map=b;}
				else if(a.equals("out")){out=b;}
				else if(a.equals("outmap")){outMap=b;}
				else if(a.equals("outfamilies")){outFamilies=b;}
				else if(a.equals("transl_table") || a.equals("table")){code=GeneticCode.forTable(Integer.parseInt(b));}
				else if(a.equals("minlen")){minLen=Parse.parseIntKMG(b);}
				else if(a.equals("cap")){cap=Parse.parseIntKMG(b);}
				else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
				else{throw new IllegalArgumentException("Unknown parameter "+arg);}
			}
			if(in==null || map==null || out==null || outMap==null || outFamilies==null){throw new IllegalArgumentException("Required: in= map= out= outmap= outfamilies=");}
			if(code==null){code=GeneticCode.forTable(11);}
			if(code.id()!=11){throw new IllegalArgumentException("This native-enumeration diagnostic only supports table 11");}
			if(minLen<3 || cap<1){throw new IllegalArgumentException("minlen and cap must be positive");}
		}
		String in, map, out, outMap, outFamilies;
		GeneticCode code;
		int minLen=80, cap=500;
		boolean overwrite;
	}

	private static final class Target{
		Target(final String id_, final String contig_, final int start1_, final int stop1_, final String strand_, final int family_){
			id=id_; contig=contig_; start1=start1_; stop1=stop1_; family=family_;
			if(strand_.length()!=1 || (strand_.charAt(0)!='+' && strand_.charAt(0)!='-')){
				throw new IllegalArgumentException("Strand must be + or - for "+id);
			}
			strand=(strand_.charAt(0)=='+' ? 0 : 1);
			if(start1<1 || stop1<start1 || family<0){throw new IllegalArgumentException("Invalid target row for "+id);}
		}
		final String id, contig;
		final int start1, stop1, strand, family;
	}

	private static final class AminoCopy{
		static byte[] reverseComplement(final byte[] bases){
			final byte[] copy=bases.clone();
			dna.AminoAcid.reverseComplementBasesInPlace(copy);
			return copy;
		}
	}

	private final Settings settings;
}
