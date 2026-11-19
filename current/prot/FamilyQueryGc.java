package prot;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.Objects;

import ddl.DDLFamilySearch.StructuralInputs;
import dna.GeneticCode;
import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.LineParser1;
import prot.cluster.FamilyTailGcSweep;
import shared.Shared;
import stream.FastaReadInputStream;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Immutable query GC in original panel order, shared by calibration and scoring.
 * The SHORT annotation and shadow factories require every query to be measured;
 * neutral amino-acid controls are selected explicitly. No inference, hashing,
 * sequence rewrite or per-query intermediate file is performed.
 * Load serially before starting scoring/read workers: native reader settings are
 * scoped and restored, but unrelated users of those globals cannot be locked here.
 * @author Brian, Fischl
 */
public final class FamilyQueryGc{

	/*--------------------------------------------------------------*/
	/*----------------            Loading           ----------------*/
	/*--------------------------------------------------------------*/

	/** Load an exact ordered panel roster from the explicit ten-field GC manifest.
	 * Validate the complete manifest before opening any panel/annotation/DNA stream.
	 * Source kinds are never inferred, and unused fields cannot hide another input.
	 * This initial interface requires structural verification and exact panel counts.
	 */
	public static synchronized FamilyQueryGc[] loadPanels(Path manifest, FamilyPanelRunner.Panel[] panels,
			InputContext inputs)throws IOException{
		require(manifest!=null && manifest.isAbsolute() && manifest.normalize().equals(manifest),
			"Query GC manifest must be an explicit normalized absolute path");
		require(panels!=null && panels.length>0 && inputs!=null,
			"Query GC manifest requires a nonempty panel roster and retained input context");
		final FamilyPanelRunner.Panel[] roster=panels.clone();
		for(FamilyPanelRunner.Panel panel:roster){
			require(panel!=null && panel.verificationMode==FamilyCalibration.VerificationMode.STRUCTURAL &&
				panel.expectedQueries>0 && panel.expectedQueries<=Integer.MAX_VALUE-8,
				"Query GC manifest loading requires structural panels with exact positive counts");
		}
		final Path manifestPath=inputs.track(manifest);
		final PanelSource[] sources=new PanelSource[roster.length];
		final FamilyQueryGc[] result=new FamilyQueryGc[roster.length];
		try(InputSettings settings=new InputSettings()){
			final ByteFile reader=ByteFile.makeByteFile1(manifestPath.toString(), false);
			final LineParser1 p=new LineParser1('\t');
			int rows=0;
			try{
				header(reader, PANEL_HEADER);
				for(byte[] row=reader.nextLine(); row!=null; row=reader.nextLine()){
					require(rows<roster.length, "Query GC manifest has extra panel rows");
					p.set(row); require(p.terms()==10, "Query GC manifest requires exactly ten fields");
					final FamilyPanelRunner.Panel panel=roster[rows];
					final Path path=absolutePath(p, 1);
					final long count=unsigned(p, 2);
					require(p.termEquals(panel.name, 0) && path.equals(panel.path) && count==panel.expectedQueries,
						"Query GC manifest panel name/path/count or original roster order differs");
					final long annotationRows=unsigned(p, 5), genomes=unsigned(p, 9);
					if(p.termEquals("annotation", 3)){
						require(annotationRows>=count && p.termEquals("-", 6) && p.termEquals("-", 7) &&
							p.termEquals("-", 8) && p.termEquals("0", 9),
							"Annotation panel requires sufficient annotation rows and unused shadow paths/count");
						sources[rows]=new PanelSource(SourceKind.ANNOTATION, absolutePath(p, 4), annotationRows, null, null, null, 0);
					}else if(p.termEquals("genomic_shadow", 3)){
						require(p.termEquals("-", 4) && p.termEquals("0", 5) && genomes>0 && genomes<=Integer.MAX_VALUE-8,
							"Genomic-shadow panel requires a positive genome count and unused annotation path/count");
						sources[rows]=new PanelSource(SourceKind.GENOMIC_SHADOW, null, 0,
							absolutePath(p, 6), absolutePath(p, 7), absolutePath(p, 8), genomes);
					}else{
						require(p.termEquals("none", 3), "Query GC manifest requires an explicit known source kind");
						require(p.termEquals("-", 4) && p.termEquals("0", 5) && p.termEquals("-", 6) && p.termEquals("-", 7) &&
							p.termEquals("-", 8) && p.termEquals("0", 9),
							"Neutral panel requires all source paths unused and both source counts zero");
						sources[rows]=new PanelSource(SourceKind.NONE, null, 0, null, null, null, 0);
					}
					rows++;
				}
				require(rows==roster.length, "Query GC manifest is missing panel rows");
			}finally{require(!reader.close(), "Query GC manifest input failed");}
			inputs.verify();
			for(int i=0; i<roster.length; i++){
				final FamilyPanelRunner.Panel panel=roster[i]; final PanelSource source=sources[i];
				if(source.kind==SourceKind.ANNOTATION){
					result[i]=annotation(panel.path, panel.expectedQueries, source.queryGc, source.annotationRows, inputs);
				}else if(source.kind==SourceKind.GENOMIC_SHADOW){
					result[i]=shadows(panel.path, panel.expectedQueries, source.sidecar, source.contigs,
						source.genomeManifest, source.genomes, inputs);
				}else{result[i]=neutral(panel.path, panel.expectedQueries, inputs);}
			}
		}
		verifyAll(result); return result;
	}

	/** Join distinct SHORT panel IDs exactly once to the accepted all-query annotation.
	 * Source indices remain annotation ordinals, independently of panel order.
	 * Counts retain the annotation producer's interval convention, including ambiguity.
	 */
	public static synchronized FamilyQueryGc annotation(Path panel, long expectedQueries, Path queryGc,
			long expectedAnnotationRows, InputContext inputs)throws IOException{
		require(expectedAnnotationRows>=expectedQueries, "Annotation population must cover the complete panel");
		final FamilyQueryGc result=new FamilyQueryGc(panel, expectedQueries, SourceKind.ANNOTATION, inputs);
		final Path panelPath=inputs.track(panel), annotationPath=inputs.track(queryGc);
		try(InputSettings settings=new InputSettings()){
			final IdIndex index=new IdIndex(result.ids.length);
			readPanel(panelPath, result, (ordinal, read)->{
				requireShort(read);
				index.add(result.ids[ordinal], ordinal);
			});
			final ByteFile reader=ByteFile.makeByteFile1(annotationPath.toString(), false);
			final LineParser1 p=new LineParser1('\t'); final ByteBuilder formatted=new ByteBuilder();
			long rows=0, matched=0;
			try{
				header(reader, QUERY_HEADER);
				for(byte[] line=reader.nextLine(); line!=null; line=reader.nextLine()){
					p.set(line);
					require(p.terms()==14 && unsigned(p, 0)==rows && rows<expectedAnnotationRows,
						"Annotation requires the declared dense original indices and fourteen fields");
					p.setBounds(1); final int i=index.get(line, p.a(), p.b());
					if(i>=0){
						if(result.sourceIndices[i]>=0){
							throw new IllegalArgumentException("Panel ID matched duplicate annotation rows: "+result.ids[i]);
						}
						require(p.termEquals("assigned", 3) || p.termEquals("unassigned", 3), "Unexpected selected annotation status");
						require(equal(p, 4, 5), "Selected SHORT annotation gene/tail GC differs");
						final long canonical=unsigned(p, 8), ambiguous=unsigned(p, 9);
						final long gc=FamilyTailGcSweep.recoverNumerator(scaled(p, 5), canonical);
						formatted.clear().append(gc/(double)canonical, 12);
						require(matches(p, 5, formatted), "Exact GC numerator must roundtrip through the original producer formatter");
						result.canonical[i]=canonical; result.ambiguous[i]=ambiguous; result.gc[i]=gc;
						result.sourceIndices[i]=rows; matched++;
					}
					rows++;
				}
			}finally{require(!reader.close(), "Annotation input failed while building query GC");}
			require(rows==expectedAnnotationRows && matched==result.ids.length, "Incomplete exact annotation/panel join");
		}
		result.verify(); return result;
	}

	/** Validate each genomic shadow's exact peptide and stop against its original DNA.
	 * Original manifest/dictionary indices and zero-based half-open intervals are retained.
	 * Only one genome's DNA is loaded at a time; peptides and join scratch are discarded.
	 */
	public static synchronized FamilyQueryGc shadows(Path panel, long expectedQueries, Path sidecar,
			Path contigs, Path genomeManifest, long expectedGenomes, InputContext inputs)throws IOException{
		final FamilyQueryGc result=new FamilyQueryGc(panel, expectedQueries, SourceKind.GENOMIC_SHADOW, inputs);
		final Path panelPath=inputs.track(panel), sidecarPath=inputs.track(sidecar);
		final Path contigPath=inputs.track(contigs), manifestPath=inputs.track(genomeManifest);
		try(InputSettings settings=new InputSettings()){
			final ArrayList<Genome> genomes=readManifest(manifestPath, expectedGenomes);
			readContigs(contigPath, genomes);
			final int[] next=new int[result.ids.length]; Arrays.fill(next, -1);
			final byte[][] peptides=new byte[result.ids.length][];
			final byte[] frames=new byte[result.ids.length];
			final IdIndex unique=new IdIndex(result.ids.length);
			final ByteFile reader=ByteFile.makeByteFile1(sidecarPath.toString(), false);
			final LineParser1 p=new LineParser1('\t');
			try{
				header(reader, SHADOW_HEADER);
				readPanel(panelPath, result, (i, read)->{
					requireShort(read); unique.add(result.ids[i], i);
					final byte[] row=reader.nextLine(); require(row!=null, "Shadow peptide has no coordinate row"); p.set(row);
					require(p.terms()==15 && p.termEquals(result.ids[i], 0),
						"Shadow sidecar ID/order or fifteen-field schema differs");
					final int gi=integer(p, 1), ci=integer(p, 2), stop=integer(p, 3), strand=integer(p, 4), frame=integer(p, 5);
					require(gi<genomes.size(), "Shadow genome index is outside the original manifest");
					final Genome genome=genomes.get(gi);
					require(ci<genome.contigs.size(), "Shadow contig index is outside its original dictionary");
					require(result.ids[i].equals("g"+gi+"_c"+ci+"_s"+stop+"_r"+strand),
						"Shadow ID does not encode its declared source coordinates");
					final int table=integer(p, 6), lo=integer(p, 10), hi=integer(p, 11);
					require((table==4 || table==11 || table==25) && integer(p, 7)==1 && integer(p, 9)==0,
						"Shadow genetic code, cohort or overlap kind differs from the bound SHORT panel");
					require(integer(p, 8)<4, "Invalid inherited genome-GC bin; this bin is not query GC");
					validateRank(p, 14);
					require(hi>lo && (long)hi-lo==3L*(read.bases.length+1),
						"Shadow interval must contain the whole peptide plus stop");
					require(strand<=1 && frame<3 && stop==(strand==0 ? hi-3 : lo), "Shadow stop orientation or frame is invalid");
					final long full=unsigned(p, 12), tail=unsigned(p, 13), span=(long)hi-lo;
					require(full>0 && full<=span && tail<=span-3 && tail<=full && full-tail<=3,
						"StopTailOverlap reports clipped overlap unions; full/tail coverage differs only within the stop");
					result.sourceIndices[i]=gi; result.contigIndices[i]=ci;
					result.starts[i]=lo; result.ends[i]=hi; result.strands[i]=(byte)strand; result.codes[i]=(byte)table;
					peptides[i]=read.bases; frames[i]=(byte)frame; next[i]=genome.first; genome.first=i;
				});
				require(reader.nextLine()==null, "Shadow sidecar has rows absent from the original panel");
			}finally{require(!reader.close(), "Shadow sidecar input failed");}
			long verified=0;
			for(Genome genome:genomes){
				if(genome.first<0){continue;}
				final Path dnaPath=inputs.track(genome.fasta);
				final HashMap<String, byte[]> dna=new HashMap<String, byte[]>();
				fasta(dnaPath, false, read->{
					require(dna.put(token(read.id), read.bases)==null, "Duplicate nucleotide contig token in "+genome.id);
				});
				final ArrayList<String> names=new ArrayList<String>(dna.keySet()); java.util.Collections.sort(names);
				require(names.equals(genome.contigs), "Actual DNA contig dictionary differs from retained source: "+genome.id);
				for(int i=genome.first; i>=0; i=next[i]){
					final byte[] bases=dna.get(genome.contigs.get(result.contigIndices[i]));
					result.gc[i]=validateDna(bases, peptides[i], result.starts[i],
						result.ends[i], result.strands[i], frames[i], result.codes[i]);
					result.canonical[i]=(long)result.ends[i]-result.starts[i];
					peptides[i]=null; verified++;
				}
			}
			require(verified==result.ids.length, "Not every original shadow query received measured nucleotide GC");
		}
		result.verify(); return result;
	}

	/** Explicit amino-acid-only panel: preserve each ordinal/ID and return neutral GC.
	 * Duplicate tokens remain distinct by ordinal; no sequence-derived GC is invented.
	 */
	public static synchronized FamilyQueryGc neutral(Path panel, long expectedQueries,
			InputContext inputs)throws IOException{
		final FamilyQueryGc result=new FamilyQueryGc(panel, expectedQueries, SourceKind.NONE, inputs);
		final Path panelPath=inputs.track(panel);
		try(InputSettings settings=new InputSettings()){readPanel(panelPath, result, (i, read)->{});}
		result.verify(); return result;
	}

	/** Allocate exclusively owned arrays; none are exposed or mutated after a factory returns. */
	private FamilyQueryGc(Path panel_, long expected, SourceKind kind_, InputContext inputs_){
		require(panel_!=null && expected>0 && expected<=Integer.MAX_VALUE-8 && kind_!=null && inputs_!=null,
			"Resident query GC requires a declared nonempty int-addressable panel and retained input context");
		final int n=(int)expected; kind=kind_; inputs=inputs_; ids=new String[n]; panel=panel_.toAbsolutePath().normalize();
		final boolean measured=kind!=SourceKind.NONE, shadow=kind==SourceKind.GENOMIC_SHADOW;
		gc=measured ? new long[n] : null; canonical=measured ? new long[n] : null;
		ambiguous=measured ? new long[n] : null; sourceIndices=measured ? new long[n] : null;
		if(measured){Arrays.fill(sourceIndices, -1);}
		contigIndices=shadow ? new int[n] : null; starts=shadow ? new int[n] : null; ends=shadow ? new int[n] : null;
		strands=shadow ? new byte[n] : null; codes=shadow ? new byte[n] : null;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Read-only worker API     ----------------*/
	/*--------------------------------------------------------------*/

	public long size(){return ids.length;}
	/** Require the original declared panel, not a different file with coincident IDs. */
	public void requirePanel(Path declared, long expectedQueries){
		require(declared!=null && panel.equals(declared.toAbsolutePath().normalize()) && expectedQueries==ids.length,
			"Query GC binding requires the exact declared panel path and query count");
	}
	/** Verify each retained context once per round, even when several panels share it. */
	public static void verifyAll(FamilyQueryGc[] panels)throws IOException{
		require(panels!=null, "Query GC verification requires explicit panel bindings");
		final IdentityHashMap<InputContext, Boolean> seen=new IdentityHashMap<InputContext, Boolean>();
		for(FamilyQueryGc panel:panels){
			require(panel!=null, "Missing query GC panel binding in verification");
			if(seen.put(panel.inputs, Boolean.TRUE)==null){panel.inputs.verify();}
		}
	}
	/** Assert the caller is scoring this panel's exact original ordinal and first token. */
	public void checkQuery(long originalOrdinal, String firstTokenId){
		final int i=index(originalOrdinal);
		if(!ids[i].equals(firstTokenId)){
			throw new IllegalArgumentException("Query GC ID/order differs at panel ordinal "+originalOrdinal);
		}
	}
	public boolean available(long ordinal){index(ordinal); return kind!=SourceKind.NONE;}
	/** Missing GC returns finite zero only together with available=false, never a measured zero. */
	public double fraction(long ordinal){
		final int i=index(ordinal); return kind==SourceKind.NONE ? 0 : gc[i]/(double)canonical[i];
	}
	public long gcBases(long ordinal){final int i=index(ordinal); return gc==null ? 0 : gc[i];}
	public long canonicalBases(long ordinal){final int i=index(ordinal); return canonical==null ? 0 : canonical[i];}
	public long ambiguousBases(long ordinal){final int i=index(ordinal); return ambiguous==null ? 0 : ambiguous[i];}
	public SourceKind sourceKind(){return kind;}
	/** Unavailable source fields return -1, governed by sourceKind; counts never use -1. */
	public long sourceIndex(long ordinal){final int i=index(ordinal); return sourceIndices==null ? -1 : sourceIndices[i];}
	public int contigIndex(long ordinal){final int i=index(ordinal); return contigIndices==null ? -1 : contigIndices[i];}
	public int start0(long ordinal){final int i=index(ordinal); return starts==null ? -1 : starts[i];}
	public int end0(long ordinal){final int i=index(ordinal); return ends==null ? -1 : ends[i];}
	public int strand(long ordinal){final int i=index(ordinal); return strands==null ? -1 : strands[i];}
	public int geneticCode(long ordinal){final int i=index(ordinal); return codes==null ? -1 : codes[i];}
	/** Check the retained declared/resolved bindings before and after downstream work. */
	public void verify()throws IOException{inputs.verify();}
	private int index(long ordinal){
		if(ordinal<0 || ordinal>=ids.length){
			throw new IllegalArgumentException("Query GC ordinal outside original panel: "+ordinal);
		}
		return (int)ordinal;
	}
	public enum SourceKind{NONE, ANNOTATION, GENOMIC_SHADOW}

	/*--------------------------------------------------------------*/
	/*----------------         Native joins          ----------------*/
	/*--------------------------------------------------------------*/

	/** Native ordered FASTA reading preserves each first token; no prefix sampling. */
	private static void readPanel(Path path, FamilyQueryGc result, PanelAction action)throws IOException{
		final int[] count={0};
		fasta(path, true, read->{
			final int i=count[0]; require(i<result.ids.length, "Panel contains more than the declared query count");
			result.ids[i]=token(read.id); action.accept(i, read); count[0]++;
		});
		require(count[0]==result.ids.length, "Panel contains fewer than the declared query count");
	}
	/** Use the native zero-worker FASTA path, with error propagation and complete closure. */
	private static void fasta(Path path, boolean amino, ReadAction action)throws IOException{
		Shared.AMINO_IN=amino;
		final FileFormat format=FileFormat.testInput(path.toString(), FileFormat.FASTA, null, true, false);
		require(format.fasta(), "Query GC requires an explicit FASTA input: "+path);
		final Streamer reader=StreamerFactory.makeStreamer(format, null, 0, true, -1, false, true, 0); reader.start();
		try{
			for(ListNum<Read> batch=reader.nextList(); batch!=null && batch.size()>0; batch=reader.nextList()){
				for(Read read:batch){
					require(read!=null && read.id!=null && read.bases!=null && read.bases.length>0 && read.mate==null,
						"Expected a nonempty unpaired FASTA record");
					action.accept(read);
				}
			}
		}finally{reader.close();}
		require(!reader.errorState(), "FASTA input failed while loading query GC: "+path);
	}
	/** Original manifest row ordinals are genome indices; retain inherited pins without hashing. */
	private static ArrayList<Genome> readManifest(Path path, long expected){
		require(expected>0 && expected<=Integer.MAX_VALUE-8, "Require a declared original genome population");
		final ArrayList<Genome> genomes=new ArrayList<Genome>(); final HashSet<String> ids=new HashSet<String>();
		final ByteFile reader=ByteFile.makeByteFile1(path.toString(), false); final LineParser1 p=new LineParser1('\t');
		try{
			header(reader, MANIFEST_HEADER);
			for(byte[] row=reader.nextLine(); row!=null; row=reader.nextLine()){
				p.set(row); require(p.terms()==9 && genomes.size()<expected, "Original genome manifest count/schema differs");
				final String id=p.parseString(0), pin=p.parseString(2); final Path fasta=Paths.get(p.parseString(1));
				require(!id.isEmpty() && ids.add(id) && pin.matches("[0-9a-f]{20}") &&
					fasta.isAbsolute() && fasta.normalize().equals(fasta),
					"Manifest requires distinct genome IDs, literal inherited identities and normalized absolute DNA paths");
				genomes.add(new Genome(id, fasta));
			}
		}finally{require(!reader.close(), "Original genome manifest input failed");}
		require(genomes.size()==expected, "Original genome manifest population differs"); return genomes;
	}
	/** Check every retained dictionary row against its original genome and dense lexical index. */
	private static void readContigs(Path path, ArrayList<Genome> genomes){
		final ByteFile reader=ByteFile.makeByteFile1(path.toString(), false); final LineParser1 p=new LineParser1('\t');
		try{
			header(reader, CONTIG_HEADER);
			for(byte[] row=reader.nextLine(); row!=null; row=reader.nextLine()){
				p.set(row); require(p.terms()==4, "Retained contig dictionary requires four fields");
				final int gi=integer(p, 0); require(gi<genomes.size(), "Dictionary genome index is outside original manifest");
				final Genome genome=genomes.get(gi); final String name=p.parseString(3);
				require(p.termEquals(genome.id, 1) && integer(p, 2)==genome.contigs.size() && !name.isEmpty(),
					"Dictionary genome identity or dense contig index differs");
				require(genome.contigs.isEmpty() || genome.contigs.get(genome.contigs.size()-1).compareTo(name)<0,
					"Dictionary contigs must be unique and lexical");
				genome.contigs.add(name);
			}
		}finally{require(!reader.close(), "Retained contig dictionary input failed");}
	}
	/** Literal translation and frame verification; ambiguous codons fail rather than becoming neutral. */
	private static long validateDna(byte[] dna, byte[] peptide, int lo, int hi, int strand, int frame, int table){
		require(dna!=null && peptide!=null && lo>=0 && hi<=dna.length && (long)hi-lo==3L*(peptide.length+1),
			"Shadow DNA interval is out of bounds");
		final int stop=strand==0 ? hi-3 : lo;
		require(frame==(strand==0 ? stop%3 : (dna.length-stop-3)%3), "Shadow frame differs from actual source DNA");
		final GeneticCode code=GeneticCode.forTable(table);
		for(int i=0; i<=peptide.length; i++){
			final int offset=strand==0 ? lo+3*i : hi-3-3*i;
			int codon=GeneticCode.codon(dna, offset);
			require(codon>=0, "Ambiguous DNA cannot supply this canonical shadow peptide");
			if(strand==1){codon=((3^(codon&3))<<4)|((3^((codon>>>2)&3))<<2)|(3^(codon>>>4));}
			final byte aa=code.aminoAcid(codon);
			require(i==peptide.length ? aa=='*' : aa==peptide[i] && aa!='*',
				"Source DNA does not reproduce the exact shadow peptide and terminal stop");
		}
		long gc=0;
		for(int i=lo; i<hi; i++){final int base=dna[i]&~32; if(base=='C' || base=='G'){gc++;}}
		return gc;
	}

	/*--------------------------------------------------------------*/
	/*----------------     Parsing and scratch       ----------------*/
	/*--------------------------------------------------------------*/

	/** Preserve explicit declared path identity; do not silently resolve relative or redundant spelling. */
	private static Path absolutePath(LineParser1 p, int field){
		final String text=p.parseString(field); final Path path=Paths.get(text);
		require(path.isAbsolute() && path.normalize().equals(path) && path.toString().equals(text),
			"Query GC source paths must be explicit normalized absolute paths");
		return path;
	}
	/** Small per-panel configuration, discarded after the validated resident tables are built. */
	private static final class PanelSource{
		PanelSource(SourceKind kind_, Path queryGc_, long annotationRows_, Path sidecar_,
				Path contigs_, Path genomeManifest_, long genomes_){
			kind=kind_;
			queryGc=queryGc_;
			annotationRows=annotationRows_;
			sidecar=sidecar_;
			contigs=contigs_;
			genomeManifest=genomeManifest_;
			genomes=genomes_;
		}
		final SourceKind kind;
		final Path queryGc, sidecar, contigs, genomeManifest;
		final long annotationRows, genomes;
	}

	/** Parse producer fractions without transient strings; the later formatter check rejects noncanonical spellings. */
	private static long scaled(LineParser1 p, int field){
		if(p.length(field)==1){
			final int digit=p.parseByte(field, 0)-'0';
			require(digit==0 || digit==1, "Measured GC requires canonical zero/one or twelve decimals");
			return digit*SCALE;
		}
		require(p.length(field)==14 && (p.parseByte(field, 0)=='0' || p.parseByte(field, 0)=='1') &&
			p.parseByte(field, 1)=='.',
			"Measured GC requires canonical zero/one or twelve decimals");
		long value=p.parseByte(field, 0)-'0';
		for(int i=2; i<14; i++){
			final int digit=p.parseByte(field, i)-'0';
			require(digit>=0 && digit<=9, "Invalid GC decimal digit");
			value=value*10+digit;
		}
		require(value<=SCALE, "Measured GC exceeds one"); return value;
	}
	/** Strict checked nonnegative decimal fields; malformed or overflowing provenance fails loudly. */
	private static long unsigned(LineParser1 p, int field){
		require(p.length(field)>0, "Missing integer field in GC provenance"); long value=0;
		for(int i=0; i<p.length(field); i++){
			final int digit=p.parseByte(field, i)-'0';
			require(digit>=0 && digit<=9, "Invalid nonnegative integer in GC provenance");
			value=Math.addExact(Math.multiplyExact(value, 10), digit);
		}
		return value;
	}
	private static int integer(LineParser1 p, int field){return Math.toIntExact(unsigned(p, field));}
	/** The producer stores its unsigned rank in a signed long; do not require nonnegative rank. */
	private static void validateRank(LineParser1 p, int field){
		require(p.length(field)>0, "Missing stored shadow rank"); final boolean negative=p.parseByte(field, 0)=='-';
		long value=0; final int start=negative ? 1 : 0; require(start<p.length(field), "Missing shadow rank digits");
		for(int i=start; i<p.length(field); i++){
			final int digit=p.parseByte(field, i)-'0'; require(digit>=0 && digit<=9, "Invalid stored shadow rank");
			value=Math.subtractExact(Math.multiplyExact(value, 10), digit);
		}
		if(!negative){Math.negateExact(value);}
	}
	private static boolean equal(LineParser1 p, int a, int b){
		if(p.length(a)!=p.length(b)){return false;}
		for(int i=0; i<p.length(a); i++){if(p.parseByte(a, i)!=p.parseByte(b, i)){return false;}}return true;
	}
	private static boolean matches(LineParser1 p, int field, ByteBuilder b){
		if(p.length(field)!=b.length){return false;}
		for(int i=0; i<b.length; i++){if(p.parseByte(field, i)!=b.array[i]){return false;}}return true;
	}
	private static void header(ByteFile reader, String expected){
		require(Arrays.equals(reader.nextLine(), expected.getBytes(StandardCharsets.US_ASCII)),
			"Unexpected GC source table header");
	}
	private static String token(String id){
		int end=0; while(end<id.length() && !Character.isWhitespace(id.charAt(end))){end++;}
		require(end>0, "FASTA header has no first token"); return id.substring(0, end);
	}
	private static void requireShort(Read read){
		require(read.bases.length>=40 && read.bases.length<=99,
			"Measured GC factory requires the bound SHORT forty-through-ninety-nine peptide population");
	}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
	private interface ReadAction{void accept(Read read)throws IOException;}
	private interface PanelAction{void accept(int ordinal, Read read)throws IOException;}

	/** Primitive ordinal slots and retained UTF-8 keys avoid allocation for unmatched annotation rows. */
	private static final class IdIndex{
		IdIndex(int n){
			require(n>0 && n<(1<<29), "Exact ID index exceeds supported resident capacity");
			int capacity=2; while(capacity<2L*n){capacity<<=1;}
			slots=new int[capacity]; keys=new byte[n][]; mask=capacity-1;
		}
		void add(String id, int ordinal){
			final byte[] key=id.getBytes(StandardCharsets.UTF_8); int slot=hash(key, 0, key.length)&mask;
			while(slots[slot]!=0){
				if(Arrays.equals(keys[slots[slot]-1], key)){
					throw new IllegalArgumentException("Duplicate measured panel query ID: "+id);
				}
				slot=(slot+1)&mask;
			}
			keys[ordinal]=key; slots[slot]=ordinal+1;
		}
		int get(byte[] bytes, int start, int end){
			int slot=hash(bytes, start, end)&mask;
			while(slots[slot]!=0){
				final int ordinal=slots[slot]-1; final byte[] key=keys[ordinal]; boolean same=key.length==end-start;
				for(int i=0; same && i<key.length; i++){same=key[i]==bytes[start+i];}
				if(same){return ordinal;}slot=(slot+1)&mask;
			}
			return -1;
		}
		/** Stable byte-range hash shared by retained keys and allocation-free annotation probes. */
		private static int hash(byte[] bytes, int start, int end){
			int value=1;
			for(int i=start; i<end; i++){value=31*value+bytes[i];}
			return value^(value>>>16);
		}
		final int[] slots;
		final byte[][] keys;
		final int mask;
	}
	/** Metadata for one original genome; the temporary linked list indexes panel-owned arrays. */
	private static final class Genome{
		Genome(String id_, Path fasta_){
			id=id_;
			fasta=fasta_;
		}
		final String id;
		final Path fasta;
		final ArrayList<String> contigs=new ArrayList<String>();
		int first=-1;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Input lifetime          ----------------*/
	/*--------------------------------------------------------------*/

	/** Caller-retained observable stability, including original symlinks and resolved targets.
	 * The shared StructuralInputs tracker accepts only regular paths, so register the
	 * actual targets there and additionally retain declared-path snapshots here.
	 * This is not a fresh digest or a guarantee against changes that preserve metadata.
	 */
	public static final class InputContext{
		public InputContext(StructuralInputs tracked_){
			require(tracked_!=null, "Query GC requires a caller-retained stability tracker");
			tracked=tracked_;
		}
		private synchronized Path track(Path path)throws IOException{
			require(path!=null, "Query GC source path is required"); final Path declared=path.toAbsolutePath().normalize();
			BoundPath bound=paths.get(declared);
			if(bound==null){bound=new BoundPath(declared); paths.put(declared, bound);}else{bound.verify();}
			tracked.track(bound.resolved); return declared;
		}
		/** Recheck both path identity and content metadata; no bulk hash or genome rescan. */
		public synchronized void verify()throws IOException{
			for(BoundPath path:paths.values()){path.verify();}tracked.verify();
		}
		private final StructuralInputs tracked;
		private final HashMap<Path, BoundPath> paths=new HashMap<Path, BoundPath>();
	}
	/** One declaration plus its regular target; retargeting to equal-size bytes is still a change. */
	private static final class BoundPath{
		BoundPath(Path path_)throws IOException{
			path=path_; declared=attributes(path);
			require(declared.isRegularFile() || declared.isSymbolicLink(),
				"Query GC input must be a regular file or inherited symlink: "+path);
			resolved=path.toRealPath(); target=attributes(resolved);
			require(target.isRegularFile() && Files.isReadable(resolved),
				"Resolved query GC input is not a readable regular file: "+path);
		}
		void verify()throws IOException{
			require(path.toRealPath().equals(resolved) && same(declared, attributes(path)),
				"Declared query GC path or symlink target changed: "+path);
			require(same(target, attributes(resolved)), "Resolved query GC input changed: "+path);
		}
		private static BasicFileAttributes attributes(Path path)throws IOException{
			return Files.readAttributes(path, BasicFileAttributes.class, LinkOption.NOFOLLOW_LINKS);
		}
		private static boolean same(BasicFileAttributes a, BasicFileAttributes b){
			return a.isRegularFile()==b.isRegularFile() && a.isSymbolicLink()==b.isSymbolicLink() && a.size()==b.size() &&
				a.lastModifiedTime().equals(b.lastModifiedTime()) && Objects.equals(a.fileKey(), b.fileKey());
		}
		final Path path, resolved;
		final BasicFileAttributes declared, target;
	}
	/** Restore reader globals and customized buffer counts even when a source fails validation. */
	private static final class InputSettings implements AutoCloseable{
		InputSettings(){
			Shared.setThreads(1); Shared.TRIM_READ_DESCRIPTION=false; FastaReadInputStream.SPLIT_READS=false;
			Read.TO_UPPER_CASE=Read.LOWER_CASE_TO_N=Read.DOT_DASH_X_TO_N=Read.FIX_HEADER=false; Read.JUNK_MODE=Read.IGNORE_JUNK;
		}
		@Override public void close(){
			Shared.AMINO_IN=amino; Shared.TRIM_READ_DESCRIPTION=trim; FastaReadInputStream.SPLIT_READS=split;
			Read.TO_UPPER_CASE=upper; Read.LOWER_CASE_TO_N=lower; Read.DOT_DASH_X_TO_N=dot; Read.FIX_HEADER=fix;
			Read.JUNK_MODE=junk; Shared.setThreads(threads); Shared.setBuffers(buffers);
		}
		final boolean amino=Shared.AMINO_IN, trim=Shared.TRIM_READ_DESCRIPTION, split=FastaReadInputStream.SPLIT_READS;
		final boolean upper=Read.TO_UPPER_CASE, lower=Read.LOWER_CASE_TO_N, dot=Read.DOT_DASH_X_TO_N, fix=Read.FIX_HEADER;
		final int junk=Read.JUNK_MODE, threads=Shared.threads(), buffers=Shared.numBuffers();
	}

	/*--------------------------------------------------------------*/
	/*----------------             State             ----------------*/
	/*--------------------------------------------------------------*/

	private final SourceKind kind;
	private final Path panel;
	private final InputContext inputs;
	private final String[] ids;
	private final long[] gc, canonical, ambiguous, sourceIndices;
	private final int[] contigIndices, starts, ends;
	private final byte[] strands, codes;
	private static final long SCALE=1000000000000L;
	private static final String PANEL_HEADER=
		"panel\tpath\tqueries\tsource_kind\tquery_gc\tannotation_rows\tsidecar\tcontigs\tgenome_manifest\tgenomes";
	private static final String QUERY_HEADER=
		"input_index\tinput_id\tstable_family_id\tstatus\tgene_gc\ttail_gc\tcontig_gc\tgenome_gc\t"+
		"gene_canonical\tgene_ambiguous\tcontig_canonical\tcontig_ambiguous\tgenome_canonical\tgenome_ambiguous";
	private static final String SHADOW_HEADER=
		"query_id\tgenome_index\tcontig_index\tstop0\tstrand\tframe\tgenetic_code\tcohort\tgc_bin\tkind\t"+
		"full_lo\tfull_hi\tfull_shadow_bp\ttail_shadow_bp\torigin_rank";
	private static final String MANIFEST_HEADER=
		"genome_id\tfasta\tfasta_sha80\tgff\tgff_sha80\tcds_rows\ttranslated_cds\tresidues\tvalid_words";
	private static final String CONTIG_HEADER="genome_index\tgenome_id\tcontig_index\tcontig_id";
}
