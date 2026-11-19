package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import dna.AminoAcid;
import dna.GeneticCode;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import gff.GffLine;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** One-genome, stop-anchored tail extraction. No calling, clustering or model loading.
 * Earlier codons outside the retained tail are deliberately not validated.
 * @author Sayu
 */
public final class StopTailCorpus {

	public static void main(String[] args) throws IOException{
		final StopTailCorpus tool=new StopTailCorpus(args);
		try{tool.process();}finally{Shared.closeStream(tool.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public StopTailCorpus(String[] args) throws IOException{
		final PreParser pp=new PreParser(args, getClass(), false); outstream=pp.outstream;
		final Parser parser=new Parser(); parser.overwrite=false;
		final HashSet<String> seen=new HashSet<String>();
		boolean tableSet=false, assumedSet=false;
		for(String arg : pp.args){
			final int eq=arg.indexOf('=');
			require(eq>0 && eq<arg.length()-1, "Expected nonempty flag=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			require(seen.add(key), "Duplicate option: "+key);
			if(key.equals("gff")){gff=value;}
			else if(key.equals("outmeta")){outmeta=value;}
			else if(key.equals("outstats")){outstats=value;}
			else if(key.equals("outaudit")){outaudit=value;}
			else if(key.equals("genome")){genome=value;}
			else if(key.equals("source")){source=value;}
			else if(key.equals("taxnamespace")){namespace=value;}
			else if(key.equals("split")){split=value;}
			else if(key.equals("taxid")){taxid=Integer.parseInt(value);}
			else if(key.equals("phylumtaxid")){phylum=Integer.parseInt(value);}
			else if(key.equals("table")){table=Integer.parseInt(value); tableSet=true;}
			else if(key.equals("tableassumed")){tableAssumed=Parse.parseBoolean(value); assumedSet=true;}
			else if(key.equals("minlen")){minlen=Integer.parseInt(value);}
			else if(key.equals("maxoriginal")){maxOriginal=Integer.parseInt(value);}
			else if(key.equals("cap")){cap=Integer.parseInt(value);}
			else if(key.equals("rejectpseudo")){rejectPseudo=Parse.parseBoolean(value);}
			else if(key.equals("reject5primepartial")){rejectFivePrime=Parse.parseBoolean(value);}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown StopTailCorpus option: "+arg);}
		}
		if(!assumedSet){tableAssumed=!tableSet;}
		tableExplicit=tableSet;
		in=parser.in1; out=parser.out1; overwrite=parser.overwrite;
		require(in!=null && gff!=null && out!=null && outmeta!=null && outstats!=null,
			"Required: in=FASTA gff=GFF out=FASTA outmeta=TSV outstats=TSV");
		require(parser.in2==null && parser.out2==null && !parser.append && parser.maxReads<0,
			"StopTailCorpus requires a whole unpaired genome; no in2/out2/append/maxreads");
		require(safe(genome) && safe(namespace), "genome and taxnamespace require nonempty ASCII letters/digits/._- only");
		require("ncbi_bacteria".equals(source) || "ncbi_archaea".equals(source) || "bioactive".equals(source)
			|| "prokcc".equals(source),
			"source must be ncbi_bacteria, ncbi_archaea, bioactive or prokcc");
		require("build".equals(split) || "calibration".equals(split) || "final".equals(split),
			"split must be build, calibration or final");
		require(taxid>=0 && phylum>=0 && table>0, "taxid/phylumtaxid must be nonnegative and table positive");
		require(minlen>=1 && cap>=minlen, "Require 1<=minlen<=cap");
		require(maxOriginal>=minlen,
			"maxoriginal must be at least minlen; it filters selected annotations, not tail truncation");
		validatePaths();
		require(Tools.testInputFiles(false, true, in, gff), "Unreadable FASTA/GFF input");
		require(Tools.testOutputFiles(overwrite, false, false, out, outmeta, outstats, outaudit),
			"Outputs must be writable and respect overwrite=f");
	}

	private void validatePaths() throws IOException{
		final String[] names=outaudit==null ? new String[]{in, gff, out, outmeta, outstats} :
			new String[]{in, gff, out, outmeta, outstats, outaudit};
		final File[] paths=new File[names.length];
		for(int i=0; i<names.length; i++){
			paths[i]=new File(names[i]).getCanonicalFile();
			// Physical files make input alias checks and three-output accounting unambiguous.
			require(!names[i].equals("-") && !names[i].startsWith("stdin") && !names[i].startsWith("stdout")
				&& !names[i].startsWith("/dev/"), "Use physical input/output files: "+names[i]);
			for(int j=0; j<i; j++){
				require(!paths[i].equals(paths[j]) && !(paths[i].exists() && paths[j].exists()
					&& Files.isSameFile(paths[i].toPath(), paths[j].toPath())),
					"Input/output paths alias: "+names[i]+" and "+names[j]);
			}
			// Keep the same fail-closed diagnostic with or without JVM assertions.
			if(i>=2){require(overwrite || !paths[i].exists(), "Outputs must be writable and respect overwrite=f: "+names[i]);}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		final HashMap<String,byte[]> bases=AnnotateGenes.readScaffolds(in);
		final TreeMap<String,Contig> contigs=new TreeMap<String,Contig>();
		int ordinal=0;
		for(String id : new TreeSet<String>(bases.keySet())){contigs.put(id, new Contig(id, ++ordinal, bases.get(id)));}
		final ArrayList<Row> rows=readGff(gff, contigs, table, tableAssumed, tableExplicit, source);
		final TreeMap<String,Row> selected=new TreeMap<String,Row>();
		for(Row row : rows){
			if(row.reason!=null){continue;}
			// Apply the explicit population policy before choosing a same-stop winner.
			if(rejectPseudo && "true".equals(attr(row.g, "pseudo"))){row.reason=Reason.PSEUDO; continue;}
			row.tail=extract(row.g, row.contig.bases, row.contig.oriented(row.g.strand), row.table, row.assumed, minlen, cap);
			row.reason=row.tail.reason;
			if(row.reason==null && rejectFivePrime && row.tail.partialFive){row.reason=Reason.FIVE_PRIME_PARTIAL;}
			if(row.reason!=null){continue;}
			row.key=genome+"|c"+row.contig.ordinal+"|"+(row.g.strand==GffLine.PLUS ? "+" : "-")
				+"|"+(row.g.strand==GffLine.PLUS ? row.g.stop : row.g.start)+"|t"+row.table;
			final Row prior=selected.get(row.key);
			if(prior==null || prefer(row, prior)){selected.put(row.key, row);}
		}
		final long[] counts=new long[Reason.values().length];
		long residues=0;
		for(Row row : rows){
			if(row.reason==null){
				final Row winner=selected.get(row.key);
				winner.aliases.add(row.nativeId);
				row.reason=row!=winner ? Reason.DUPLICATE_ALIAS :
					row.tail.originalLength>maxOriginal ? Reason.TOO_LONG : Reason.EMITTED;
				if(row.reason==Reason.EMITTED){residues=Math.addExact(residues, row.tail.peptide.length);}
			}
			counts[row.reason.ordinal()]++;
		}
		// Preserve the same longest-eligible stop winner before restricting the output band.
		// Filtering before dedup would promote short aliases of otherwise valid long genes.
		if(maxOriginal<Integer.MAX_VALUE){
			final java.util.Iterator<Row> winners=selected.values().iterator();
			while(winners.hasNext()){
				if(winners.next().reason==Reason.TOO_LONG){winners.remove();}
			}
		}
		long accounted=0; for(long count : counts){accounted=Math.addExact(accounted, count);}
		require(accounted==rows.size() && counts[Reason.EMITTED.ordinal()]==selected.size(),
			"Every CDS needs exactly one accounting outcome");
		write(selected, rows, counts, residues, contigs.size());
		outstream.println("STOP_TAIL_CORPUS_PASS cds="+rows.size()+" emitted="+selected.size()
			+" duplicate_alias="+counts[Reason.DUPLICATE_ALIAS.ordinal()]+" residues="+residues);
	}

	static GenomeView parseGenomeForControls(String fasta, String gff, String source,
			int fallbackTable, boolean tableAssumed, boolean tableExplicit){
		validateControlInputs(fasta, gff, source, fallbackTable);
		final HashMap<String,byte[]> bases=AnnotateGenes.readScaffolds(fasta);
		final TreeMap<String,Contig> contigs=new TreeMap<String,Contig>();
		int ordinal=0;
		for(String id : new TreeSet<String>(bases.keySet())){contigs.put(id, new Contig(id, ++ordinal, bases.get(id)));}
		final ArrayList<Row> rows=readGff(gff, contigs, fallbackTable, tableAssumed, tableExplicit, source);
		return new GenomeView(contigs, rows);
	}

	private static void validateControlInputs(String fasta, String gff, String source, int fallbackTable){
		require(fasta!=null && gff!=null, "Control parser requires FASTA and GFF inputs");
		require(!fasta.equals(gff), "Control parser FASTA/GFF paths must be distinct");
		require("ncbi_bacteria".equals(source) || "ncbi_archaea".equals(source) || "bioactive".equals(source)
			|| "prokcc".equals(source),
			"Control parser source must be ncbi_bacteria, ncbi_archaea, bioactive or prokcc");
		require(fallbackTable>0, "Control parser fallback translation table must be positive");
		require(Tools.testInputFiles(false, true, fasta, gff), "Control parser inputs must be readable");
	}

	static int annotatedStop1(Row row){
		require(row.g.strand==GffLine.PLUS || row.g.strand==GffLine.MINUS,
			"Annotated stop masks require a known strand at GFF line "+row.lineNumber);
		final int endpoint=row.g.strand==GffLine.PLUS ? row.g.stop : row.g.start;
		// Validated excluded circular children can have either endpoint in virtual space.
		return ((endpoint-1)%row.contig.bases.length)+1;
	}

	static ArrayList<Row> readGff(String gff, TreeMap<String,Contig> contigs,
			int table, boolean tableAssumed, boolean tableExplicit, String source){
		readTopology(gff, contigs);
		final ArrayList<Row> rows=new ArrayList<Row>();
		final HashMap<String,Group> groups=new HashMap<String,Group>();
		final ArrayList<GffLine> parentFeatures=new ArrayList<GffLine>();
		final ByteFile input=ByteFile.makeByteFile(gff, false);
		final LineParser1 fields=new LineParser1('\t');
		final boolean oldSource=GffLine.paseSource;
		GffLine.paseSource=true;
		int fallback=table, lineNumber=0;
		boolean assumed=tableAssumed, features=false, directive=false, shifted=false;
		final boolean predicted="bioactive".equals(source) || "prokcc".equals(source);
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lineNumber++;
				if(line.length==0){continue;}
				if(line[0]=='#'){
					final String text=new String(line, StandardCharsets.UTF_8);
					require(!text.equals("##FASTA"), "Embedded GFF FASTA is unsupported; supply in= separately");
					if(text.startsWith("##transl_table ") || text.startsWith("##transl_table\t")
							|| text.startsWith("##genetic-code custom")){
						require(!features && !directive, "Translation-code directive must occur once before features");
						fallback=text.startsWith("##genetic-code custom") ? -1 : Integer.parseInt(text.substring(14).trim());
						require(fallback==-1 || fallback>0, "Invalid file translation table");
						assumed=predicted && tableAssumed; directive=true;
					}
					continue;
				}
				features=true; fields.set(line);
				require(fields.terms()==9, "Expected nine GFF columns at line "+lineNumber);
				for(int i=0; i<9; i++){require(!fields.termEquals("", i), "Empty GFF field at line "+lineNumber);}
				final long start=fields.parseLong(3), stop=fields.parseLong(4);
				require(start>=1 && stop>=start && stop<=Integer.MAX_VALUE, "Invalid GFF coordinates at line "+lineNumber);
				require(fields.termEquals('+', 6) || fields.termEquals('-', 6) || fields.termEquals('.', 6)
					|| fields.termEquals('?', 6),
					"Invalid GFF strand at line "+lineNumber);
				require(fields.termEquals('.', 7) || fields.termEquals('0', 7) || fields.termEquals('1', 7)
					|| fields.termEquals('2', 7),
					"Invalid GFF phase at line "+lineNumber);
				final GffLine g=new GffLine(line);
				final Contig contig=contigs.get(AnnotationGff.decode(g.seqid));
				require(contig!=null, "GFF interval outside/missing FASTA contig at line "+lineNumber);
				final int length=contig.bases.length;
				final boolean wrap=g.stop>length, cds="CDS".equals(g.type);
				// NCBI circular-origin intervals legitimately extend past the linear FASTA.
				// V1 excludes these CDSs rather than clipping or inventing wrapped translation.
				// Shifted NCBI CDS children are provisional until parent/companion validation below.
				final boolean shiftedChild=cds && ("ncbi_bacteria".equals(source) || "ncbi_archaea".equals(source))
					&& (g.strand==GffLine.PLUS || g.strand==GffLine.MINUS);
				require(!wrap || contig.topology==1 && (g.start<=length || shiftedChild) && (long)g.stop<=2L*length
					&& (long)g.stop-g.start+1<=length, "GFF interval outside/missing FASTA contig at line "+lineNumber);
				if(!cds){
					if("gene".equals(g.type) || "mRNA".equals(g.type)){parentFeatures.add(g);}
					continue;
				}
				shifted|=g.start>length;
				g.partial();
				partialEnds(g); // Malformed ranges fail even when a multipart group is later excluded.
				final String declared=attr(g, "transl_table"), id=attr(g, "ID"), parent=attr(g, "Parent");
				final int code=declared==null ? fallback : Integer.parseInt(declared);
				require(code==-1 || code>0, "Invalid CDS translation table at line "+lineNumber);
				final boolean codeAssumed=declared==null ? assumed : predicted && tableAssumed;
				final String codeSource=declared!=null ? "feature" : directive ? "file" : tableExplicit ? "genome" : "default";
				final Row row=new Row(g, contig, code, codeAssumed, codeSource,
					id==null ? "." : AnnotationGff.decode(id), lineNumber);
				if(wrap){row.reason=Reason.CIRCULAR_ORIGIN;}
				rows.add(row);
				addGroup(groups, id==null ? null : "ID:"+AnnotationGff.decode(id), row);
				if(parent!=null){
					for(int from=0; from<parent.length(); ){
						int end=parent.indexOf(',', from); if(end<0){end=parent.length();}
						require(end>from, "Empty CDS Parent at line "+lineNumber);
						addGroup(groups, "Parent:"+AnnotationGff.decode(parent.substring(from, end)), row);
						from=end+1; require(from!=parent.length(), "Trailing empty CDS Parent at line "+lineNumber);
					}
				}
			}
		}finally{GffLine.paseSource=oldSource; require(!input.close(), "Failed reading GFF: "+gff);}
		if(shifted){
			final HashMap<String,GffLine> parents=new HashMap<String,GffLine>();
			// Ordinary inputs never needed parent attributes; do not impose new parsing failures on them.
			for(GffLine parent : parentFeatures){
				final String id=attr(parent, "ID");
				if(id!=null){
					final String key=AnnotationGff.decode(id);
					// Null permanently marks an ambiguous parent, including a non-wrapping duplicate.
					if(parents.containsKey(key)){parents.put(key, null);}else{parents.put(key, parent);}
				}
			}
			for(Row row : rows){if(row.g.start>row.contig.bases.length){validateShiftedChild(row, groups, parents);}}
		}
		for(Group group : groups.values()){
			if(group.compound){for(Row row : group.rows){if(row.reason==null){row.reason=Reason.COMPOUND_CDS;}}}
		}
		return rows;
	}

	/** Recognize the observed Pandoraea multipart representation, not arbitrary start>L GFF.
	 * NCBI's documented origin convention does not explicitly promise shifted child starts.
	 * Require a real single-turn parent and crossing same-ID companion; never translate these rows. */
	private static void validateShiftedChild(Row row, HashMap<String,Group> groups, HashMap<String,GffLine> parents){
		final int length=row.contig.bases.length;
		assert(row.g.start>length && row.reason==Reason.CIRCULAR_ORIGIN) :
			"Only provisional excluded circular children reach validation";
		final String message="GFF interval outside/missing FASTA contig at line "+row.lineNumber
			+": unsupported shifted circular CDS child";
		final String parent=attr(row.g, "Parent");
		require(!row.nativeId.isEmpty() && !row.nativeId.equals(".") && parent!=null && !parent.isEmpty()
			&& !parent.equals(".") && parent.indexOf(',')<0, message+"; need CDS ID and one Parent");
		final String parentId=AnnotationGff.decode(parent);
		final GffLine anchor=parents.get(parentId);
		require(anchor!=null && row.contig.id.equals(AnnotationGff.decode(anchor.seqid)) && anchor.strand==row.g.strand
			&& anchor.start<=length && anchor.stop>length && (long)anchor.stop<=2L*length
			&& (long)anchor.stop-anchor.start+1<=length, message+"; need unique real origin-crossing gene/mRNA parent");
		final Group group=groups.get("ID:"+row.nativeId);
		require(group!=null && group.rows.size()>1, message+"; need same-ID companion");
		boolean crossing=false;
		for(Row part : group.rows){
			final String linked=attr(part.g, "Parent");
			require(part.contig==row.contig && part.g.strand==row.g.strand && linked!=null && linked.indexOf(',')<0
				&& parentId.equals(AnnotationGff.decode(linked)) && part.g.start>=anchor.start && part.g.stop<=anchor.stop,
				message+"; all same-ID segments must share parent/contig/strand and lie within parent");
			crossing|=part.g.start<=length && part.g.stop>length;
		}
		// Containment in a <=L parent also bounds the complete group to one turn.
		require(crossing, message+"; need same-ID origin-crossing companion");
	}

	/** Read explicit topology independently of feature order; never infer it from coordinates. */
	private static void readTopology(String gff, TreeMap<String,Contig> contigs){
		final ByteFile input=ByteFile.makeByteFile(gff, false);
		final LineParser1 fields=new LineParser1('\t');
		int lineNumber=0;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lineNumber++;
				if(line.length==0 || line[0]=='#'){continue;}
				fields.set(line);
				require(fields.terms()==9, "Expected nine GFF columns at line "+lineNumber);
				if(!fields.termEquals("region", 2)){continue;}
				final String declaration=CdsTranslator.attribute(fields.parseString(8), "Is_circular");
				if(declaration==null){continue;}
				require(declaration.equals("true") || declaration.equals("false"),
					"Invalid Is_circular declaration at line "+lineNumber);
				final Contig contig=contigs.get(AnnotationGff.decode(fields.parseString(0)));
				require(contig!=null, "GFF interval outside/missing FASTA contig at line "+lineNumber);
				require(fields.parseLong(3)==1 && fields.parseLong(4)==contig.bases.length,
					"Topology region must cover the whole FASTA contig at line "+lineNumber);
				final int topology=declaration.equals("true") ? 1 : 2;
				require(contig.topology==0 || contig.topology==topology,
					"Conflicting Is_circular declarations at line "+lineNumber);
				if(contig.topology==0){
					contig.topology=topology; contig.topologyLine=lineNumber;
					contig.topologyDeclaration=new String(line, StandardCharsets.UTF_8);
				}
			}
		}finally{require(!input.close(), "Failed reading GFF topology: "+gff);}
	}

	private static void addGroup(HashMap<String,Group> groups, String key, Row row){
		if(key==null){return;}
		Group group=groups.get(key);
		if(group==null){group=new Group(); groups.put(key, group);}
		if(!group.rows.isEmpty()){
			final Row first=group.rows.get(0);
			final boolean same=first.contig==row.contig && first.g.start==row.g.start
				&& first.g.stop==row.g.stop && first.g.strand==row.g.strand;
			final boolean alternate=key.startsWith("Parent:") && first.contig==row.contig
				&& first.g.strand==row.g.strand && (row.g.strand==GffLine.PLUS && first.g.stop==row.g.stop
				|| row.g.strand==GffLine.MINUS && first.g.start==row.g.start)
				&& !first.nativeId.equals(".") && !row.nativeId.equals(".") && !first.nativeId.equals(row.nativeId);
			group.compound|=!same && !alternate;
		}
		group.rows.add(row);
	}

	private static boolean prefer(Row a, Row b){
		assert(a.key.equals(b.key)) : "Only alternative annotations of the same biological stop/code may compete";
		if(a.tail.originalLength!=b.tail.originalLength){return a.tail.originalLength>b.tail.originalLength;}
		int order=a.nativeId.compareTo(b.nativeId);
		if(order!=0){return order<0;}
		order=a.g.attributes.compareTo(b.g.attributes);
		if(order!=0){return order<0;}
		if(a.g.start!=b.g.start){return a.g.start<b.g.start;}
		if(a.g.stop!=b.g.stop){return a.g.stop<b.g.stop;}
		if(a.g.phase!=b.g.phase){return a.g.phase<b.g.phase;}
		return a.g.source.compareTo(b.g.source)<0;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Tail semantics        ----------------*/
	/*--------------------------------------------------------------*/

	/** Public literal-fixture entry; the command caches one reverse complement per contig. */
	public static Tail extract(GffLine g, byte[] bases, int fallback, boolean assumed, int minlen, int cap){
		require(g!=null && bases!=null, "CDS and FASTA bases are required");
		final byte[] oriented=g.strand==GffLine.MINUS ? bases.clone() : bases;
		if(g.strand==GffLine.MINUS){AminoAcid.reverseComplementBasesInPlace(oriented);}
		final String declared=attr(g, "transl_table");
		return extract(g, bases, oriented, declared==null ? fallback : Integer.parseInt(declared),
			declared==null && assumed, minlen, cap);
	}

	static Tail extract(GffLine g, byte[] bases, byte[] oriented, int table, boolean assumed, int minlen, int cap){
		require(g.start>=1 && g.stop>=g.start && g.stop<=bases.length, "CDS coordinates must lie inside FASTA");
		require(minlen>=1 && cap>=minlen, "Require 1<=minlen<=cap");
		final Tail result=new Tail(table, assumed);
		if(g.strand!=GffLine.PLUS && g.strand!=GffLine.MINUS){return result.reject(Reason.UNKNOWN_STRAND);}
		if(g.phase<0){return result.reject(Reason.UNKNOWN_PHASE);}
		require(g.phase<=2, "CDS phase must be 0,1,2 or unknown");
		final int partial=partialEnds(g);
		if(partial<0){return result.reject(Reason.AMBIGUOUS_PARTIAL);}
		result.partialFive=(partial&1)!=0; result.partialThree=(partial&2)!=0;
		if(result.partialThree){return result.reject(Reason.THREE_PRIME_PARTIAL);}
		if(attr(g, "transl_except")!=null || attr(g, "exception")!=null){return result.reject(Reason.TRANSLATION_EXCEPTION);}
		if(table!=4 && table!=11 && table!=25){return result.reject(Reason.UNSUPPORTED_TRANSLATION_TABLE);}
		final int length=g.stop-g.start+1-g.phase;
		if(length<3 || length%3!=0){return result.reject(Reason.INCOMPLETE_TERMINAL_CODON);}
		result.originalLength=length/3-1;
		final int start=(g.strand==GffLine.PLUS ? g.start-1 : bases.length-g.stop)+g.phase;
		final int stop=g.strand==GffLine.PLUS ? g.stop-1 : bases.length-g.start;
		final GeneticCode code=GeneticCode.forTable(table);
		if(!code.isStop(GeneticCode.codon(oriented, stop-2))){return result.reject(Reason.MISSING_TERMINAL_STOP);}
		if(result.originalLength<minlen){return result.reject(Reason.TOO_SHORT);}
		result.peptide=ProteinSupport.translateTail(oriented, start, stop, code, cap);
		for(byte amino : result.peptide){
			if(amino=='*'){return result.reject(Reason.INTERNAL_STOP);}
			if(amino=='X'){return result.reject(Reason.AMBIGUOUS_RESIDUE);}
		}
		assert(result.peptide.length==Math.min(cap, result.originalLength)) :
			"Stop triplet must be omitted and cap applied at the upstream edge only";
		return result;
	}

	/** Return biological five/three bits, or -1 when the declared sides conflict/are unknown. */
	static int partialEnds(GffLine g){
		g.partial();
		final String partial=attr(g, "partial"), low=attr(g, "start_range"), high=attr(g, "end_range");
		final int leftRange=rangeSide(low, g.start), rightRange=rangeSide(high, g.stop);
		if(leftRange<0 || rightRange<0){return -1;}
		boolean left=leftRange==1, right=rightRange==1;
		if("true".equals(partial)){
			if(!left && !right){return -1;}
		}else if(partial!=null){
			final boolean pl=partial.length()==2 && partial.charAt(0)=='1';
			final boolean pr=partial.length()==2 && partial.charAt(1)=='1';
			if((low!=null && left!=pl) || (high!=null && right!=pr)){return -1;}
			left=pl; right=pr;
		}
		return g.strand==GffLine.MINUS ? (right ? 1 : 0)|(left ? 2 : 0) : (left ? 1 : 0)|(right ? 2 : 0);
	}

	private static int rangeSide(String value, int coordinate){
		if(value==null){return 0;}
		final int comma=value.indexOf(',');
		require(comma>0 && comma<value.length()-1 && value.indexOf(',', comma+1)<0, "Malformed partial range: "+value);
		final String a=value.substring(0, comma), b=value.substring(comma+1);
		final long low=a.equals(".") ? -1 : Long.parseLong(a), high=b.equals(".") ? -1 : Long.parseLong(b);
		require((a.equals(".") || low>0) && (b.equals(".") || high>0), "Invalid partial range: "+value);
		if(low==-1 && high==coordinate || high==-1 && low==coordinate){return 1;}
		return low==coordinate && high==coordinate ? 0 : -1;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(TreeMap<String,Row> selected, ArrayList<Row> rows, long[] counts, long residues, int contigs){
		ByteStreamWriter fasta=null, meta=null, stats=null, audit=null;
		boolean failed=false;
		final ByteBuilder line=new ByteBuilder();
		try{
			fasta=ByteStreamWriter.makeBSW(out, overwrite, false, false);
			meta=ByteStreamWriter.makeBSW(outmeta, overwrite, false, false);
			stats=ByteStreamWriter.makeBSW(outstats, overwrite, false, false);
			if(outaudit!=null){audit=ByteStreamWriter.makeBSW(outaudit, overwrite, false, false);}
			meta.println("tail_id\tgenome_id\tcontig\tcontig_ordinal\tstrand\tstop_coordinate"
				+"\tcds_start\tcds_end\tphase\tnative_id\tnative_aliases\tsource\tgff_source"
				+"\ttaxid\tphylum_taxid\ttaxonomy_namespace\ttable\tcode_source\ttable_assumed"
				+"\toriginal_length\ttail_length\tpartial_5prime\tpartial_3prime\tsplit");
			for(Row row : selected.values()){
				final Tail tail=row.tail; final GffLine g=row.g;
				line.clear().append('>').append(row.key).append(" taxid=").append(taxid).append(" phylum_taxid=").append(phylum)
					.append(" taxnamespace=").append(namespace).append(" source=").append(source).append(" split=").append(split);
				fasta.println(line); fasta.println(tail.peptide);
				line.clear().append(row.key).tab().append(genome).tab(); escaped(line, row.contig.id);
				line.tab().append(row.contig.ordinal).tab().append(g.strand==GffLine.PLUS ? '+' : '-').tab()
					.append(g.strand==GffLine.PLUS ? g.stop : g.start).tab().append(g.start).tab()
					.append(g.stop).tab().append(g.phase).tab();
				escaped(line, row.nativeId); line.tab();
				boolean first=true;
				for(String alias : row.aliases){if(!first){line.append(',');} escaped(line, alias); first=false;}
				line.tab().append(source).tab(); escaped(line, g.source);
				line.tab().append(taxid).tab().append(phylum).tab().append(namespace).tab()
					.append(tail.table).tab().append(row.codeSource).tab()
					.append(tail.assumed ? 1 : 0).tab().append(tail.originalLength).tab().append(tail.peptide.length).tab()
					.append(tail.partialFive ? 1 : 0).tab().append(tail.partialThree ? 1 : 0).tab().append(split);
				meta.println(line);
			}
			if(audit!=null){
				audit.println("genome_id\tgff_line\tcontig\tstart\tend\tstrand\tphase\tnative_id"
					+"\tstatus\treason\ttail_id\ttable\tcode_source\ttable_assumed\tcontig_length"
					+"\ttopology_gff_line\ttopology_declaration");
				for(Row row : rows){
					line.clear().append(genome).tab().append(row.lineNumber).tab(); escaped(line, row.contig.id);
					line.tab().append(row.g.start).tab().append(row.g.stop).tab()
						.append(row.g.strand==GffLine.PLUS ? '+' : row.g.strand==GffLine.MINUS ? '-' : '.')
						.tab().append(row.g.phase).tab(); escaped(line, row.nativeId);
					line.tab().append(row.reason==Reason.EMITTED ? "emitted" :
						row.reason==Reason.DUPLICATE_ALIAS ? "duplicate_alias" : "excluded")
						.tab().append(row.reason.name().toLowerCase(Locale.ROOT)).tab().append(row.key==null ? "." : row.key)
						.tab().append(row.table).tab().append(row.codeSource).tab().append(row.assumed ? 1 : 0);
					line.tab().append(row.contig.bases.length).tab().append(row.contig.topologyLine).tab();
					escaped(line, row.contig.topologyDeclaration);
					audit.println(line);
				}
			}
			stats.println("schema\tgenome_id\tmetric\tcount");
			metric(stats, "contigs", contigs);
			metric(stats, "cds_total", rows.size());
			metric(stats, "tail_residues", residues);
			metric(stats, "excluded", rows.size()-selected.size()-counts[Reason.DUPLICATE_ALIAS.ordinal()]);
			for(Reason reason : Reason.values()){
				if(reason!=Reason.TOO_LONG || maxOriginal<Integer.MAX_VALUE){
					metric(stats, reason.name().toLowerCase(Locale.ROOT), counts[reason.ordinal()]);
				}
			}
			final TreeMap<Integer,long[]> codes=new TreeMap<Integer,long[]>();
			for(Row row : rows){
				long[] count=codes.get(row.table); if(count==null){count=new long[2]; codes.put(row.table, count);}
				count[row.assumed ? 1 : 0]++;
			}
			for(Map.Entry<Integer,long[]> entry : codes.entrySet()){
				metric(stats, "table_"+entry.getKey()+"_declared", entry.getValue()[0]);
				metric(stats, "table_"+entry.getKey()+"_assumed", entry.getValue()[1]);
			}
		}finally{
			if(fasta!=null){failed|=fasta.poisonAndWait();}
			if(meta!=null){failed|=meta.poisonAndWait();}
			if(stats!=null){failed|=stats.poisonAndWait();}
			if(audit!=null){failed|=audit.poisonAndWait();}
			require(!failed, "Failed writing tail FASTA, sidecar or accounting");
		}
	}

	private void metric(ByteStreamWriter out, String key, long value){
		out.println(new ByteBuilder().append("stop_tail_v1").tab().append(genome).tab()
			.append(key).tab().append(value));
	}
	private static void escaped(ByteBuilder out, String value){AnnotationGff.escaped(out, value==null ? "." : value);}
	private static String attr(GffLine g, String name){return CdsTranslator.attribute(g.attributes, name);}
	private static boolean safe(String value){return value!=null && value.matches("[A-Za-z0-9_.-]+");}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}

	public enum Reason {
		EMITTED, DUPLICATE_ALIAS, COMPOUND_CDS, TRANSLATION_EXCEPTION, UNSUPPORTED_TRANSLATION_TABLE,
		UNKNOWN_STRAND, UNKNOWN_PHASE, AMBIGUOUS_PARTIAL, THREE_PRIME_PARTIAL, INCOMPLETE_TERMINAL_CODON,
		MISSING_TERMINAL_STOP, TOO_SHORT, INTERNAL_STOP, AMBIGUOUS_RESIDUE, CIRCULAR_ORIGIN, TOO_LONG,
		PSEUDO, FIVE_PRIME_PARTIAL
	}
	public static final class Tail {
		Tail(int table_, boolean assumed_){table=table_; assumed=assumed_;}
		Tail reject(Reason reason_){reason=reason_; peptide=null; return this;}
		public final int table;
		public final boolean assumed;
		public boolean partialFive, partialThree;
		public int originalLength;
		public byte[] peptide;
		public Reason reason;
	}
	static final class GenomeView {
		GenomeView(TreeMap<String,Contig> contigs_, ArrayList<Row> rows_){contigs=contigs_; rows=rows_;}
		final TreeMap<String,Contig> contigs;
		final ArrayList<Row> rows;
	}
	static final class Contig {
		Contig(String id_, int ordinal_, byte[] bases_){id=id_; ordinal=ordinal_; bases=bases_;}
		byte[] oriented(int strand){
			if(strand!=GffLine.MINUS){return bases;}
			if(reverse==null){reverse=bases.clone(); AminoAcid.reverseComplementBasesInPlace(reverse);}
			return reverse;
		}
		final String id;
		final int ordinal;
		final byte[] bases;
		byte[] reverse;
		int topology, topologyLine;
		String topologyDeclaration;
	}
	static final class Row {
		Row(GffLine g_, Contig contig_, int table_, boolean assumed_, String codeSource_,
				String nativeId_, int lineNumber_){
			g=g_; contig=contig_; table=table_; assumed=assumed_; codeSource=codeSource_;
			nativeId=nativeId_; lineNumber=lineNumber_;
		}
		final GffLine g;
		final Contig contig;
		final int table;
		final boolean assumed;
		final String nativeId;
		final String codeSource;
		final int lineNumber;
		final TreeSet<String> aliases=new TreeSet<String>();
		String key;
		Tail tail;
		Reason reason;
	}
	private static final class Group {final ArrayList<Row> rows=new ArrayList<Row>(); boolean compound;}
	private final PrintStream outstream;
	private String in, gff, out, outmeta, outstats, outaudit, genome, source, namespace, split;
	private int taxid=-1, phylum=-1, table=11, minlen=100, cap=200;
	private int maxOriginal=Integer.MAX_VALUE;
	private boolean tableAssumed=true, tableExplicit, overwrite;
	private boolean rejectPseudo, rejectFivePrime;
}
