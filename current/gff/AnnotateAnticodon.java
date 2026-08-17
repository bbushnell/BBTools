package gff;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;

/**
 * Annotates tRNA features in a GFF with the anticodon TRIPLET sequence.
 *
 * NCBI GFFs record a tRNA's anticodon by genomic POSITION (anticodon=(pos:X..Y) or
 * anticodon=(pos:complement(X..Y))) and the amino acid by name (product=tRNA-Xxx), but NOT the
 * triplet as text. This tool reads the 3 bases at that position from the paired genome FASTA,
 * reverse-complements them when the tRNA is on the minus strand, converts DNA-&gt;RNA (T-&gt;U), and
 * appends Note=tRNA-Xxx(YYY) to the feature's attributes. That Note= form is the first-priority
 * branch of TrnaConsensusBuilder.parseAnticodon and is carried into FASTA headers by CutGff, so an
 * annotated GFF flows through cutgff -&gt; TrnaConsensusBuilder with zero downstream code changes.
 *
 * tRNAscan-SE GFF3 (-j) instead states the amino acid as isotype=Xxx and the anticodon DIRECTLY as a
 * triplet (anticodon=ACG, already in tRNA 5'-&gt;3' orientation). Those forms are also recognized: the
 * triplet is used as-is (uppercased, T-&gt;U) with no genome lookup or strand revcomp, so tRNAscan output
 * is normalized to the same Note=tRNA-Xxx(YYY) key and groups consistently with the position-annotated
 * corpus. product=/anticodon=(pos:...) take priority, so mixed or NCBI-style GFFs are unaffected.
 *
 * All input lines pass through UNCHANGED except tRNA lines, which gain the Note= attribute. tRNA
 * features already carrying Note=tRNA-Xxx(YYY) are left untouched (idempotent). A tRNA with no
 * anticodon position is labeled Note=tRNA-Xxx(UNK) (not skipped), matching extract_trnas.py's ac='UNK'
 * fallback so all such tRNAs cluster under one UNK group. This is the Java reimplementation of
 * scripts/extract_trnas.py's extraction logic (lines 104-115); it matches that behavior, including
 * skipping (not annotating) only a feature whose seqid is absent from the FASTA.
 *
 * Usage: annotateanticodon.sh in=genome.fna.gz gff=genome.gff.gz out=annotated.gff.gz
 *   Batch: in=a.fna.gz,b.fna.gz gff=a.gff.gz,b.gff.gz outdir=staging/  (gff= inferred if omitted)
 *
 * @author Brian Bushnell, Noire, Amber
 * @date August 13, 2026
 */
public class AnnotateAnticodon {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Command-line entry point.
	 * @param args Command line arguments */
	public static void main(String[] args){
		Timer t=new Timer();
		AnnotateAnticodon x=new AnnotateAnticodon(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructs the tool from command-line arguments.
	 * @param args Command line arguments
	 */
	public AnnotateAnticodon(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		Shared.TRIM_READ_DESCRIPTION=Shared.TRIM_RNAME=true;//keep only the first header token as Read.id
		Read.TO_UPPER_CASE=true;

		parse(args);

		fnaList=Tools.fixExtension(fnaList);
		gffList=Tools.fixExtension(gffList);
		inferGff();
		checkFileExistence();
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses command-line arguments.
	 * @param args Array of command line arguments */
	private void parse(String[] args){
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("in") || a.equals("fna") || a.equals("infna") || a.equals("ref")){
				assert(b!=null) : "in= requires a value";
				Tools.addFiles(b, fnaList);
			}else if(a.equals("gff") || a.equals("ingff")){
				assert(b!=null) : "gff= requires a value";
				Tools.addFiles(b, gffList);
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("outdir")){
				outdir=b;
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(arg.indexOf('=')<0 && new File(arg).exists()){
				fnaList.add(arg);
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		overwrite=parser.overwrite;
		append=parser.append;
	}

	/** Infers a .gff.gz path from each .fna path when gff= was not supplied. */
	private void inferGff(){
		if(!gffList.isEmpty()){
			assert(gffList.size()==fnaList.size()) : "fna/gff count mismatch: "+fnaList.size()+" vs "+gffList.size();
			return;
		}
		for(String s : fnaList){
			String prefix=ReadWrite.stripExtension(s);
			String gff=prefix+".gff";
			if(!new File(gff).exists()){
				String gz=gff+".gz";
				assert(new File(gz).exists()) : "Can't find gff for "+s+" (tried "+gff+" and "+gz+")";
				gff=gz;
			}
			gffList.add(gff);
		}
	}

	/** Validates inputs and output target. */
	private void checkFileExistence(){
		if(fnaList.isEmpty()){throw new RuntimeException("Error - at least one input (in=) is required.");}
		assert(gffList.size()==fnaList.size()) : "fna/gff count mismatch: "+fnaList.size()+" vs "+gffList.size();
		ArrayList<String> in=new ArrayList<String>();
		in.addAll(fnaList); in.addAll(gffList);
		if(!Tools.testInputFiles(false, true, in.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read some input files.\n");
		}
		if(outdir==null && out==null){throw new RuntimeException("Error - out= or outdir= is required.");}
		if(outdir==null && fnaList.size()>1){
			throw new RuntimeException("Error - multiple inputs require outdir= (a directory), not out=.");
		}
		if(outdir!=null){
			File d=new File(outdir);
			if(!d.exists()){d.mkdirs();}
			assert(d.isDirectory()) : "outdir is not a directory: "+outdir;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Annotates every input file pair.
	 * @param t Timer for reporting */
	void process(Timer t){
		for(int i=0; i<fnaList.size(); i++){
			String fna=fnaList.get(i), gff=gffList.get(i);
			String outPath=(outdir!=null ? outdir+"/"+new File(gff).getName() : out);
			annotateFile(fna, gff, outPath);
		}
		t.stop();
		outstream.println("Files:             \t"+fnaList.size());
		outstream.println("tRNA lines:        \t"+trnaLines);
		outstream.println("Annotated:         \t"+annotated);
		outstream.println("Already annotated: \t"+alreadyAnnotated);
		outstream.println("No anticodon pos:  \t"+noAnticodonPos);
		outstream.println("Scaffold not found:\t"+scaffoldNotFound);
		outstream.println("Bad anticodon pos: \t"+badAnticodonPos);
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, bytesProcessed, 8));
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Annotates one FASTA/GFF pair, writing the modified GFF to outPath.
	 * @param fna Genome FASTA path
	 * @param gff Input GFF path
	 * @param outPath Output (annotated) GFF path
	 */
	private void annotateFile(String fna, String gff, String outPath){
		final HashMap<String, byte[]> genome=loadFasta(fna);

		final FileFormat ffgff=FileFormat.testInput(gff, FileFormat.GFF, null, true, true);
		final FileFormat ffout=FileFormat.testOutput(outPath, FileFormat.GFF, null, true, overwrite, append, false);
		final ByteFile bf=ByteFile.makeByteFile(ffgff);
		final ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();

		final ByteBuilder bb=new ByteBuilder();
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			linesProcessed++;
			bytesProcessed+=(line.length+1);
			if(line.length==0){bsw.println(line); continue;}
			if(line[0]=='#'){bsw.println(line); continue;}

			//Cheap type gate: only fully parse tRNA lines
			if(!isTrnaLine(line)){bsw.println(line); continue;}

			trnaLines++;
			String annotatedLine=annotateTrnaLine(new String(line), genome, bb);
			bsw.println(annotatedLine);
		}
		errorState|=bf.close();
		bsw.poisonAndWait();
		if(verbose){outstream.println("Wrote "+outPath);}
	}

	/** True if the 3rd tab-separated field of a GFF line equals "tRNA". */
	private static boolean isTrnaLine(byte[] line){
		int tabs=0, start=-1, end=-1;
		for(int i=0; i<line.length; i++){
			if(line[i]=='\t'){
				tabs++;
				if(tabs==2){start=i+1;}
				else if(tabs==3){end=i; break;}
			}
		}
		if(start<0 || end<0 || end-start!=4){return false;}//"tRNA" is 4 chars
		return line[start]=='t' && line[start+1]=='R' && line[start+2]=='N' && line[start+3]=='A';
	}

	/**
	 * Returns the tRNA line with Note=tRNA-Xxx(YYY) appended, or the line unchanged when it cannot be
	 * (or need not be) annotated. Counters record each outcome.
	 * @param line The full tab-separated tRNA GFF line
	 * @param genome seqid -&gt; uppercase bases (with _tid_-stripped aliases)
	 * @param bb Scratch buffer
	 * @return The (possibly annotated) line
	 */
	private String annotateTrnaLine(String line, HashMap<String, byte[]> genome, ByteBuilder bb){
		String[] f=line.split("\t", -1);
		if(f.length<9){return line;}
		final String seqid=f[0], strand=f[6], attrs=f[8];

		if(NOTE_PAT.matcher(attrs).find()){alreadyAnnotated++; return line;}//idempotent

		final Matcher aaM=PRODUCT_PAT.matcher(attrs);
		final String amino;
		if(aaM.find()){amino=aaM.group(1);}
		else{final Matcher isoM=ISOTYPE_PAT.matcher(attrs); amino=(isoM.find() ? isoM.group(1) : "Unk");}//tRNAscan-SE isotype= form

		//Determine the anticodon triplet. A tRNA with no genomic anticodon position (or an out-of-bounds
		//one) is labeled (UNK) rather than skipped, matching extract_trnas.py's ac='UNK' fallback, so the
		//consensus builder groups them together under one UNK key instead of scattering them by amino acid.
		String triplet;
		final Matcher posM=ANTICODON_POS_PAT.matcher(attrs);
		if(posM.find()){//NCBI form: anticodon given by genomic position -> read bases + revcomp-by-strand + T->U
			byte[] bases=genome.get(seqid);
			if(bases==null){bases=genome.get(stripTid(seqid));}
			if(bases==null){scaffoldNotFound++; return line;}//can't resolve scaffold; skip (post-repair this is 0)
			triplet=tripletRNA(bases, Integer.parseInt(posM.group(1)), Integer.parseInt(posM.group(2)), "-".equals(strand));
			if(triplet==null){triplet="UNK"; badAnticodonPos++;}
			else{annotated++;}
		}else{
			//tRNAscan-SE form: anticodon stated directly as a triplet (already tRNA-oriented) -> use as-is, T->U.
			final Matcher seqM=ANTICODON_SEQ_PAT.matcher(attrs);
			if(seqM.find()){triplet=toRNA(seqM.group(1)); annotated++;}
			else{triplet="UNK"; noAnticodonPos++;}
		}

		bb.clear();
		for(int i=0; i<8; i++){bb.append(f[i]).tab();}
		if(".".equals(attrs) || attrs.isEmpty()){bb.append("Note=tRNA-").append(amino).append('(').append(triplet).append(')');}
		else{//attrs already present -> add a ';' separator only if it lacks a trailing one (tRNAscan attrs end in ';')
			bb.append(attrs);
			if(attrs.charAt(attrs.length()-1)!=';'){bb.append(';');}
			bb.append("Note=tRNA-").append(amino).append('(').append(triplet).append(')');
		}
		return bb.toString();
	}

	/**
	 * Extracts the anticodon triplet as RNA. Reads the 1-based inclusive genomic range acStart..acStop
	 * forward, reverse-complements it when on the minus strand, uppercases, and maps T-&gt;U.
	 * @param bases Uppercase scaffold bases (0-based)
	 * @param acStart 1-based inclusive anticodon start
	 * @param acStop 1-based inclusive anticodon stop
	 * @param minus True if the tRNA is on the minus strand
	 * @return The RNA triplet, or null if the coordinates are out of bounds
	 */
	private static String tripletRNA(byte[] bases, int acStart, int acStop, boolean minus){
		if(acStart<1 || acStop<acStart || acStop>bases.length){return null;}
		final int len=acStop-acStart+1;
		final byte[] tri=new byte[len];
		if(minus){//reverse then complement (revcomp keyed on strand, matching extract_trnas.py)
			for(int i=0; i<len; i++){tri[i]=AminoAcid.baseToComplementExtended[bases[acStop-1-i]];}
		}else{
			for(int i=0; i<len; i++){tri[i]=bases[acStart-1+i];}
		}
		for(int i=0; i<len; i++){
			byte c=tri[i];
			if(c>='a' && c<='z'){c-=32;}
			if(c=='T'){c='U';}
			tri[i]=c;
		}
		return new String(tri);
	}

	/** Uppercases a DNA anticodon triplet and maps T-&gt;U, WITHOUT reverse-complementing (tRNAscan-SE
	 * already reports the anticodon in tRNA 5'-&gt;3' orientation). @param s DNA triplet. @return RNA triplet. */
	private static String toRNA(String s){
		final byte[] b=s.getBytes();
		for(int i=0; i<b.length; i++){
			byte c=b[i];
			if(c>='a' && c<='z'){c-=32;}
			if(c=='T'){c='U';}
			b[i]=c;
		}
		return new String(b);
	}

	/**
	 * Loads a genome FASTA into a seqid-&gt;bases map, registering both the full header key and its
	 * _tid_&lt;taxid&gt;-stripped form so bare-accession GFF seqids resolve (extract_trnas.py lines 23-27).
	 * @param fna FASTA path
	 * @return Map from sequence id to uppercase bases
	 */
	private HashMap<String, byte[]> loadFasta(String fna){
		ArrayList<Read> reads=ReadInputStream.toReads(fna, FileFormat.FA, -1);
		HashMap<String, byte[]> map=new HashMap<String, byte[]>(reads.size()*4);
		for(Read r : reads){
			map.put(r.id, r.bases);
			String stripped=stripTid(r.id);
			if(!stripped.equals(r.id) && !map.containsKey(stripped)){map.put(stripped, r.bases);}
		}
		return map;
	}

	/** Removes a trailing _tid_&lt;digits&gt; suffix, or returns the input unchanged. */
	private static String stripTid(String s){
		int idx=s.lastIndexOf("_tid_");
		if(idx<0){return s;}
		final int digitStart=idx+5;
		if(digitStart>=s.length()){return s;}
		for(int i=digitStart; i<s.length(); i++){
			if(!Character.isDigit(s.charAt(i))){return s;}
		}
		return s.substring(0, idx);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<String> fnaList=new ArrayList<String>();
	private ArrayList<String> gffList=new ArrayList<String>();
	private String out=null;
	private String outdir=null;

	private long linesProcessed=0, bytesProcessed=0;
	private long trnaLines=0, annotated=0, alreadyAnnotated=0;
	private long noAnticodonPos=0, scaffoldNotFound=0, badAnticodonPos=0;

	/** Already-annotated tRNA: Note=tRNA-Xxx(YYY). */
	private static final Pattern NOTE_PAT=Pattern.compile("Note=tRNA-\\w+\\(\\w+\\)");
	/** Reference amino acid: product=tRNA-Xxx. */
	private static final Pattern PRODUCT_PAT=Pattern.compile("product=tRNA-(\\w+)");
	/** Anticodon genomic position: anticodon=(pos:X..Y) or anticodon=(pos:complement(X..Y)). */
	private static final Pattern ANTICODON_POS_PAT=Pattern.compile("anticodon=\\(pos:(?:complement\\()?(\\d+)\\.\\.(\\d+)");
	/** Reference amino acid, tRNAscan-SE form: isotype=Xxx (fallback used when product=tRNA-Xxx is absent). */
	private static final Pattern ISOTYPE_PAT=Pattern.compile("isotype=(\\w+)");
	/** Anticodon stated directly as a triplet, tRNAscan-SE form: anticodon=ACG (NOT the (pos:...) form). */
	private static final Pattern ANTICODON_SEQ_PAT=Pattern.compile("anticodon=([ACGTUacgtu]{3})(?![ACGTUacgtu])");

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;

}
