package prot;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashSet;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import structures.ByteBuilder;

/**
 * Strict, data-free contract for the shared synthetic-bin selection manifest.
 *
 * <p>One row describes one stable bin identity.  The row records the exact
 * native and foreign contig instances selected, rather than only a seed or a
 * requested fraction.  Consequently a vectorizer can replay a row after the
 * sampler changes.  The manifest deliberately permits the same organism in
 * both train and validation; validation is a random row sample, not an
 * organism holdout.</p>
 *
 * <p>Schema v1/v2 columns are identical:
 * {@code bin_id, split, target_tid, contaminant_tids, comp_requested,
 * cont_requested, spike_class, breakpoints, native_contigs, foreign_contigs,
 * model_rows}.  Contig lists use {@code tid|contig_id|instance_id} tokens separated by
 * semicolons; repeated source spans remain distinct through the instance ID.
 * The first pipe separates the tid and the last pipe separates the instance
 * ID.  Contaminant tids use semicolon-separated integers.  A dash means empty.
 *
 * <p><b>v1 (legacy):</b> raw contig IDs must not contain the list delimiters (tab, comma,
 * CR, LF) -- rejected, not escaped. v1 files are read forever exactly as before, byte-for-byte;
 * this class never applies any decoding to a v1 field, so a literal {@code %}-looking substring
 * in a v1 identifier stays completely literal.
 *
 * <p><b>v2 (2026-09-09, Ady/UMP45, Yoimiya's format-repair assignment):</b> real contig/shred
 * identifiers can contain commas and semicolons (UMP45's real-corpus gate: 547,044 comma-bearing,
 * 626 semicolon-bearing ids), which v1 cannot represent at all. v2 writers {@link IdentifierCodec#encode}
 * every contig_id/instance_id component before embedding it in a token (so the structural
 * delimiters -- tab, {@code ;}, {@code |} -- never collide with real identifier bytes, which are
 * always escaped wherever they'd otherwise collide); v2 readers must {@link IdentifierCodec#decode}
 * each component back to its true value via {@link #decodeSelections(String,int)} -- never via the
 * legacy single-argument {@link #decodeSelections(String)}, which always means "treat as v1" and is
 * kept solely so callers that have not yet been updated keep their exact current (raw) behavior
 * rather than silently mis-decoding v2 data. The schema-version header is the ONLY encoding gate
 * (no separate flag); an unrecognized version fails loud at load time.</p>
 */
public final class MagQCBinManifest {

	public static final int SCHEMA_VERSION=1;
	/** v2 (Ady/UMP45, 2026-09-09): identifiers are {@link IdentifierCodec}-encoded; see class javadoc. */
	public static final int SCHEMA_VERSION_ENCODED=2;
	private static final String HEADER="#bin_id\tsplit\ttarget_tid\tcontaminant_tids\tcomp_requested\t"
		+"cont_requested\tspike_class\tbreakpoints\tnative_contigs\tforeign_contigs\tmodel_rows";

	private MagQCBinManifest(){ }

	public static void main(String[] args){
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){selftest(); return;}
		if(args.length!=1){throw new RuntimeException("Usage: java -ea prot.MagQCBinManifest selftest|<manifest.tsv>");}
		final Manifest m=load(args[0]);
		System.err.println("MagQCBinManifest: validated "+m.size()+" bins; manifest_sha256="+sha256File(args[0]));
	}

	public static final class Bin {
		public final String id, split, spikeClass, breakpoints, nativeContigs, foreignContigs, modelRows;
		public final int targetTid;
		public final String contaminantTids;
		public final double compFraction, contFraction;
		/** 1 (legacy, raw) or 2 ({@link IdentifierCodec}-encoded) -- describes THIS Bin's OWN
		 *  nativeContigs/foreignContigs strings. Carried on the object itself (UMP45's proposal,
		 *  2026-09-09) rather than threaded through every call chain, because a Bin can legitimately
		 *  escape its parent Manifest's scope (e.g. {@link MagQCReferenceCdsIndex#labelsFor}, which
		 *  only ever sees a Bin) -- a caller with no Manifest in scope still has the right version. */
		public final int schemaVersion;

		/** Legacy convenience: always constructs a v1 (raw) Bin, matching every pre-2026-09-09 call
		 *  site's existing behavior unchanged. A caller that has (or is producing) v2-encoded
		 *  components must use {@link #Bin(String,String,int,String,double,double,String,String,String,String,String,int)}. */
		public Bin(String id, String split, int targetTid, String contaminantTids,
				double compFraction, double contFraction, String spikeClass, String breakpoints,
				String nativeContigs, String foreignContigs, String modelRows){
			this(id,split,targetTid,contaminantTids,compFraction,contFraction,spikeClass,breakpoints,
				nativeContigs,foreignContigs,modelRows,SCHEMA_VERSION);
		}

		/** Explicit-version constructor. {@code schemaVersion} must be {@link #SCHEMA_VERSION} or
		 *  {@link #SCHEMA_VERSION_ENCODED}; this class does not encode/decode anything itself --
		 *  the caller is responsible for having already {@link IdentifierCodec#encode}d
		 *  nativeContigs/foreignContigs components before calling this constructor when passing
		 *  {@link #SCHEMA_VERSION_ENCODED} (structural validation below is unaffected either way:
		 *  it only parses pipe/tid POSITIONS, never contig_id content). */
		public Bin(String id, String split, int targetTid, String contaminantTids,
				double compFraction, double contFraction, String spikeClass, String breakpoints,
				String nativeContigs, String foreignContigs, String modelRows, int schemaVersion){
			if(schemaVersion!=SCHEMA_VERSION && schemaVersion!=SCHEMA_VERSION_ENCODED){
				throw new RuntimeException("Unrecognized schema version for Bin "+id+": "+schemaVersion);
			}
			this.schemaVersion=schemaVersion;
			this.id=requireToken(id,"bin_id"); this.split=requireSplit(split); this.targetTid=targetTid;
			this.contaminantTids=normalizeList(contaminantTids,"contaminant_tids",false);
			this.compFraction=checkFraction(compFraction,"comp_requested");
			this.contFraction=checkFraction(contFraction,"cont_requested");
			this.spikeClass=requireToken(spikeClass,"spike_class");
			this.breakpoints=normalizeList(breakpoints,"breakpoints",true);
			this.nativeContigs=normalizeList(nativeContigs,"native_contigs",true);
			this.foreignContigs=normalizeList(foreignContigs,"foreign_contigs",true);
			this.modelRows=normalizeList(modelRows,"model_rows",true);
			validateSelections();
		}

		private void validateSelections(){
			final HashSet<String> nativeSet=new HashSet<String>();
			for(String token : splitList(nativeContigs)){
				if(token.length()==0){continue;}
				final int bar=token.indexOf('|');
				if(bar<=0 || bar==token.length()-1){throw new RuntimeException("Malformed native contig token '"+token+"' in "+id);}
				final int bar2=token.lastIndexOf('|');
				if(bar2<=bar+1 || bar2==token.length()-1){throw new RuntimeException("Contig token must be tid|contig_id|instance_id: '"+token+"' in "+id);}
				final int tid=parseTid(token.substring(0,bar));
				if(tid!=targetTid){throw new RuntimeException("Native contig tid "+tid+" != target_tid "+targetTid+" in "+id);}
				if(schemaVersion==SCHEMA_VERSION_ENCODED){
					checkEncodedComponent(token.substring(bar+1,bar2), "native", id);
					checkEncodedComponent(token.substring(bar2+1), "native", id);
				}
				if(!nativeSet.add(token)){throw new RuntimeException("Duplicate native contig '"+token+"' in "+id);}
			}
			final HashSet<String> all=new HashSet<String>(nativeSet);
			for(String token : splitList(foreignContigs)){
				if(token.length()==0){continue;}
				final int bar=token.indexOf('|');
				if(bar<=0 || bar==token.length()-1){throw new RuntimeException("Malformed foreign contig token '"+token+"' in "+id);}
				final int bar2=token.lastIndexOf('|');
				if(bar2<=bar+1 || bar2==token.length()-1){throw new RuntimeException("Contig token must be tid|contig_id|instance_id: '"+token+"' in "+id);}
				final int tid=parseTid(token.substring(0,bar));
				if(tid==targetTid){throw new RuntimeException("Foreign contig has target tid "+tid+" in "+id);}
				if(schemaVersion==SCHEMA_VERSION_ENCODED){
					checkEncodedComponent(token.substring(bar+1,bar2), "foreign", id);
					checkEncodedComponent(token.substring(bar2+1), "foreign", id);
				}
				if(!all.add(token)){throw new RuntimeException("Contig selected twice ('"+token+"') in "+id);}
			}
			final HashSet<String> tids=new HashSet<String>();
			for(String token : splitList(contaminantTids)){
				if(token.length()==0){continue;}
				final int tid=parseTid(token);
				if(tid==targetTid){throw new RuntimeException("Contaminant tid equals target tid in "+id);}
				tids.add(token);
			}
			for(String token : splitList(foreignContigs)){
				final String tid=token.substring(0, token.indexOf('|'));
				if(!tids.contains(tid)){throw new RuntimeException("Foreign contig tid "+tid+" absent from contaminant_tids in "+id);}
			}
		}

		/** v2 boundary validation (Yoimiya source-review finding, 2026-09-09): a malformed
		 *  percent escape must fail HERE, at construction/load time, not later at first actual
		 *  {@link IdentifierCodec#decode} inside some downstream consumer -- crash-loud at the
		 *  input boundary, not deep in the pipeline. */
		private void checkEncodedComponent(String component, String role, String binId){
			if(!IdentifierCodec.isEncoded(component)){
				throw new RuntimeException("Malformed "+IdentifierCodec.ENCODING+" escape in "+role+" contig component '"
					+component+"' of bin "+binId);
			}
		}
	}

	public static final class Manifest {
		public final String sourcePath, sourceSha256, secondarySourcePath, secondarySourceSha256;
		public final ArrayList<Bin> bins;
		/** 1 (legacy, raw) or 2 ({@link IdentifierCodec}-encoded) -- the ONLY encoding gate; thread
		 *  this through to every {@link #decodeSelections(String,int)} call a consumer makes on
		 *  this manifest's bins, rather than assuming a version. */
		public final int schemaVersion;
		/** Zero denotes a historical manifest without a declared cross-model split policy. */
		public final int splitModulus;
		Manifest(String sourcePath, String sourceSha256, String secondarySourcePath, String secondarySourceSha256,
				ArrayList<Bin> bins, int schemaVersion){
			this(sourcePath,sourceSha256,secondarySourcePath,secondarySourceSha256,bins,schemaVersion,0);
		}
		Manifest(String sourcePath,String sourceSha256,String secondarySourcePath,String secondarySourceSha256,
				ArrayList<Bin> bins,int schemaVersion,int splitModulus){
			this.sourcePath=sourcePath; this.sourceSha256=sourceSha256;
			this.secondarySourcePath=secondarySourcePath; this.secondarySourceSha256=secondarySourceSha256; this.bins=bins;
			this.schemaVersion=schemaVersion;
			this.splitModulus=splitModulus;
		}
		public int size(){return bins.size();}
	}

	/**
	 * Selects one model's rows without changing any bin definition or split.
	 * The existing model_rows column holds semicolon-separated memberships of
	 * the form {@code model:train:ordinal} or {@code model:val:ordinal}. Ordinals
	 * are zero-based and contiguous within each model/split. Legacy
	 * {@code train:ordinal}/{@code val:ordinal} tokens may coexist and do not
	 * select a named model. A bin cannot cross its manifest train/validation
	 * split. Returned entries are the original immutable Bin objects.
	 *
	 * @throws IllegalArgumentException for duplicate, missing, or malformed row ordinals
	 */
	public static ArrayList<Bin> selectModelRows(Manifest manifest,String split,String model){
		if(manifest==null || !("train".equals(split) || "val".equals(split))){
			throw new IllegalArgumentException("Model selection requires a manifest and train/val split");
		}
		validateModelName(model);
		final ArrayList<Bin> chosen=new ArrayList<Bin>();
		final structures.IntList ordinals=new structures.IntList();
		for(Bin bin : manifest.bins){
			if(!bin.split.equals(split)){continue;}
			final int ordinal=modelOrdinal(bin,model);
			if(ordinal>=0){chosen.add(bin); ordinals.add(ordinal);}
		}
		final Bin[] ordered=new Bin[chosen.size()];
		for(int i=0; i<chosen.size(); i++){
			final int ordinal=ordinals.get(i);
			if(ordinal>=ordered.length || ordered[ordinal]!=null){
				throw new IllegalArgumentException("Noncontiguous or duplicate row "+ordinal+
					" for model "+model+" split "+split+" (selected "+ordered.length+" bins)");
			}
			ordered[ordinal]=chosen.get(i);
		}
		chosen.clear();
		for(Bin bin : ordered){
			assert(bin!=null) : "Unique ordinals in [0,selected count) must fill every model row";
			chosen.add(bin);
		}
		return chosen;
	}

	/** Model IDs remain literal safe tokens inside the colon/semicolon membership grammar. */
	public static void validateModelName(String model){
		if(model==null || model.isEmpty()){
			throw new IllegalArgumentException("binmodel= must name one model");
		}
		for(int i=0; i<model.length(); i++){
			final char ch=model.charAt(i);
			if(!((ch>='A' && ch<='Z') || (ch>='a' && ch<='z') || (ch>='0' && ch<='9') ||
					ch=='_' || ch=='-' || ch=='.')){
				throw new IllegalArgumentException("Invalid character in model ID: "+model);
			}
		}
	}

	/** Parses membership ranges directly, avoiding split arrays/temporary strings per bin. */
	static int modelOrdinal(Bin bin,String model){
		final String text=bin.modelRows;
		if(text.equals("-")){return -1;}
		int selected=-1;
		for(int start=0; start<text.length(); ){
			final int separator=text.indexOf(';',start),end=(separator<0 ? text.length() : separator);
			final int colon=text.indexOf(':',start);
			if(colon<=start || colon>=end){throw new IllegalArgumentException("Malformed model_rows in "+bin.id);}
			final int second=text.indexOf(':',colon+1);
			final boolean named=(second>=0 && second<end);
			final int splitStart=(named ? colon+1 : start),splitEnd=(named ? second : colon);
			if(splitEnd-splitStart!=bin.split.length() || !text.regionMatches(splitStart,bin.split,0,bin.split.length())){
				throw new IllegalArgumentException("Model membership changes the train/validation split of "+bin.id);
			}
			final int numberStart=splitEnd+1;
			if(numberStart==end){throw new IllegalArgumentException("Missing model row ordinal in "+bin.id);}
			int ordinal=0;
			for(int j=numberStart; j<end; j++){
				final int digit=text.charAt(j)-'0';
				if(digit<0 || digit>9 || ordinal>(Integer.MAX_VALUE-digit)/10){
					throw new IllegalArgumentException("Invalid model row ordinal in "+bin.id);
				}
				ordinal=ordinal*10+digit;
			}
			if(named && colon-start==model.length() && text.regionMatches(start,model,0,model.length())){
				if(selected>=0){throw new IllegalArgumentException("Duplicate membership for model "+model+" in "+bin.id);}
				selected=ordinal;
			}
			start=end+1;
		}
		return selected;
	}

	/** Decoded explicit contig selection; the instance ID is never inferred from the contig ID. */
	public static final class Selection {
		public final int tid;
		public final String contigId, instanceId;
		Selection(int tid_, String contigId_, String instanceId_){tid=tid_; contigId=contigId_; instanceId=instanceId_;}
	}

	/** Decodes a native/foreign selection list for replay consumers, ALWAYS as v1 (legacy, raw --
	 *  no {@link IdentifierCodec} decoding applied). Kept exactly as it always behaved so every
	 *  existing caller's behavior is unchanged; a caller that knows the manifest's real
	 *  {@link Manifest#schemaVersion} must call {@link #decodeSelections(String,int)} instead, or
	 *  a v2 field's components stay wrongly still-encoded. */
	public static ArrayList<Selection> decodeSelections(String encoded){
		return decodeSelections(encoded, SCHEMA_VERSION);
	}

	/** Bin-aware convenience (UMP45's proposal, 2026-09-09): decodes {@code bin}'s native ({@code
	 *  foreign=false}) or foreign ({@code foreign=true}) selection list using the VERSION CARRIED
	 *  ON THE BIN ITSELF, never a version threaded separately through the call chain -- correct
	 *  even when the caller has no parent {@link Manifest} in scope. Prefer this over the raw
	 *  {@link #decodeSelections(String,int)} whenever a {@code Bin} is available. */
	public static ArrayList<Selection> decodeSelections(Bin bin, boolean foreign){
		return decodeSelections(foreign ? bin.foreignContigs : bin.nativeContigs, bin.schemaVersion);
	}

	/** Decodes a native/foreign selection list for replay consumers, honoring the manifest's real
	 *  schema version: {@code schemaVersion==1} splits and returns each component exactly as
	 *  written (today's behavior, unchanged -- a literal {@code %}-looking substring stays
	 *  literal); {@code schemaVersion==2} additionally {@link IdentifierCodec#decode}s the
	 *  contig_id and instance_id of every token back to their true values, AFTER the raw
	 *  tab/semicolon/pipe structural split -- never the other order, since a decoded component
	 *  can legitimately contain a literal {@code ;} or {@code |} that must never be re-split.
	 *  Any other value throws (unrecognized schema version). */
	public static ArrayList<Selection> decodeSelections(String encoded, int schemaVersion){
		if(schemaVersion!=SCHEMA_VERSION && schemaVersion!=SCHEMA_VERSION_ENCODED){
			throw new RuntimeException("Unrecognized schema version for decodeSelections: "+schemaVersion);
		}
		final ArrayList<Selection> out=new ArrayList<Selection>();
		if(encoded==null || encoded.equals("-")){return out;}
		for(String token : encoded.split(";",-1)){
			final int bar=token.indexOf('|'), bar2=token.lastIndexOf('|');
			if(bar<=0 || bar2<=bar+1 || bar2==token.length()-1){
				throw new RuntimeException("Malformed selection token: "+token);
			}
			final int tid=parseTid(token.substring(0,bar));
			String contigId=token.substring(bar+1,bar2), instanceId=token.substring(bar2+1);
			if(schemaVersion==SCHEMA_VERSION_ENCODED){
				contigId=IdentifierCodec.decode(contigId); instanceId=IdentifierCodec.decode(instanceId);
			}
			out.add(new Selection(tid, contigId, instanceId));
		}
		return out;
	}

	/** Loads a complete manifest, retaining the historical global duplicate-ID check. */
	public static Manifest load(String file){
		final HashSet<String> ids=new HashSet<String>();
		try(Reader reader=new Reader(file)){
			final Manifest manifest=reader.header();
			for(Bin bin=reader.next(); bin!=null; bin=reader.next()){
				if(!ids.add(bin.id)){throw new RuntimeException("Duplicate bin_id "+bin.id+" in "+file);}
				manifest.bins.add(bin);
			}
			return manifest;
		}
	}

	/**
	 * Streaming manifest reader: retains only headers and one pending byte line.
	 * Each returned Bin receives the same row/selection validation as load().
	 * Global ID uniqueness and per-model ordinal coverage belong to the consumer;
	 * enforcing them here would defeat bounded-memory merge passes. Use load()
	 * when a complete in-memory manifest with global ID validation is needed.
	 * Always use try-with-resources: an exception parsing next() does not close
	 * the reader automatically, and an intentionally consumed prefix has no EOF.
	 */
	public static final class Reader implements AutoCloseable {

		/** Reads and validates all headers before exposing source bindings to the caller. */
		public Reader(String file){
			this.file=file;
			input=ByteFile.makeByteFile(file,true);
			try{pending=readRow();}
			catch(RuntimeException e){closeAfterFailure(e); throw e;}
			catch(Error e){closeAfterFailure(e); throw e;}
		}

		/** Returns an independent empty manifest carrying the validated source and policy headers. */
		public Manifest header(){
			assert(sawSchema && sawSource && sawColumns) : "Reader construction must finish header validation";
			return new Manifest(sourcePath,sourceSha256,secondaryPath,secondaryHash,
				new ArrayList<Bin>(),schemaVersion,splitModulus);
		}

		/** Returns the next validated bin, or null at EOF; no returned bin is retained internally. */
		public Bin next(){
			if(eof){return null;}
			if(closed){throw new IllegalStateException("Manifest reader is closed: "+file);}
			if(pending==null){pending=readRow();}
			if(pending==null){return null;}
			parser.set(pending); pending=null;
			final Bin bin=parseBin(parser,schemaVersion,file);
			rowsRead++;
			return bin;
		}

		/** Number of rows successfully parsed, excluding headers and comments. */
		public long rowsRead(){return rowsRead;}

		/** Stops the underlying ByteFile even when a merge intentionally consumes only a prefix. */
		@Override public void close(){
			if(!closed){
				closed=true; pending=null;
				if(input.close()){throw new RuntimeException("I/O error reading "+file);}
			}
		}

		/** Constructor failures must not leave a ByteFile producer running without a consumer. */
		private void closeAfterFailure(Throwable failure){
			try{close();}catch(RuntimeException closeFailure){failure.addSuppressed(closeFailure);}
		}

		/** Consumes headers/comments and returns one data line without allocating field Strings. */
		private byte[] readRow(){
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0){continue;}
				parser.set(line);
				if(line[0]=='#'){
					readHeader(line);
					continue;
				}
				if(!sawColumns && parser.termEquals("bin_id",0)){
					acceptColumns(new String(line,java.nio.charset.StandardCharsets.UTF_8));
					continue;
				}
				requireHeaders(); dataStarted=true;
				return line;
			}
			requireHeaders(); eof=true; close();
			return null;
		}

		/** Known metadata may not change after data starts; unrelated comment lines remain legal. */
		private void readHeader(byte[] line){
			final boolean known=parser.termEquals("#schema_version",0) || parser.termEquals("#source",0) ||
				parser.termEquals("#source_secondary",0) || parser.termEquals("#realization_split",0) ||
				parser.termEquals("#columns",0) || parser.termEquals("#bin_id",0);
			if(!known){return;}
			if(dataStarted){throw new RuntimeException("Late manifest header in "+file);}
			if(parser.termEquals("#schema_version",0)){
				if(sawSchema || parser.terms()!=2){throw new RuntimeException("Duplicate or malformed schema header in "+file);}
				schemaVersion=Integer.parseInt(parser.parseString(1));
				if(schemaVersion!=SCHEMA_VERSION && schemaVersion!=SCHEMA_VERSION_ENCODED){
					throw new RuntimeException("Unrecognized schema_version "+schemaVersion+" in "+file);
				}
				sawSchema=true;
			}else if(parser.termEquals("#source",0) || parser.termEquals("#source_secondary",0)){
				final boolean secondary=parser.termEquals("#source_secondary",0);
				if(parser.terms()!=5 || !parser.termEquals("sha256",3) || (secondary ? sawSecondary : sawSource)){
					throw new RuntimeException("Duplicate or malformed source header in "+file);
				}
				if(secondary){secondaryPath=parser.parseString(1); secondaryHash=parser.parseString(4); sawSecondary=true;}
				else{sourcePath=parser.parseString(1); sourceSha256=parser.parseString(4); sawSource=true;}
			}else if(parser.termEquals("#realization_split",0)){
				if(splitModulus!=0 || parser.terms()!=4 || !parser.termEquals(MagQCBinSplit.VERSION,1) ||
						!parser.termStartsWith("modulus=",2) || !parser.termEquals("validation_bucket=0",3)){
					throw new RuntimeException("Unsupported realization split policy in "+file);
				}
				splitModulus=Integer.parseInt(parser.parseString(2).substring(8));
				MagQCBinSplit.validate(splitModulus);
			}else{
				final String text=new String(line,java.nio.charset.StandardCharsets.UTF_8);
				acceptColumns(parser.termEquals("#columns",0) ? text.substring(9) : text.substring(1));
			}
		}

		/** Accepts the documented header spellings once, with the exact eleven-column layout. */
		private void acceptColumns(String columns){
			if(sawColumns || !columns.equals(HEADER.substring(1))){throw new RuntimeException("Duplicate or invalid columns header in "+file);}
			sawColumns=true;
		}

		/** Empty manifests still require complete headers and a nonempty source binding. */
		private void requireHeaders(){
			if(!sawSchema || !sawSource || !sawColumns){throw new RuntimeException("Incomplete manifest headers in "+file);}
			if(sourceSha256.isEmpty()){throw new RuntimeException("Manifest is missing #source provenance in "+file);}
		}

		private final String file;
		private final ByteFile input;
		private final LineParser1 parser=new LineParser1((byte)'\t');
		private byte[] pending;
		private String sourcePath="-",sourceSha256="-",secondaryPath="-",secondaryHash="-";
		private int schemaVersion=-1,splitModulus=0;
		private long rowsRead;
		private boolean sawSchema,sawSource,sawSecondary,sawColumns,dataStarted,eof,closed;
	}

	/** Writes a v1 (legacy) manifest and a hash sidecar from bytes read back from disk. Kept
	 *  exactly as it always behaved -- every existing caller's Bin objects are assumed to hold
	 *  raw (unencoded) identifiers, matching the v1 contract, so this never touches
	 *  {@link IdentifierCodec}. A caller that pre-encoded its Bin objects' identifiers must call
	 *  {@link #write(String,String,String,ArrayList,int)} with {@link #SCHEMA_VERSION_ENCODED}
	 *  instead, or the file's declared version would not match its actual bytes. */
	public static String write(String file, String sourcePath, String sourceSha256, ArrayList<Bin> bins){
		return write(file, sourcePath, sourceSha256, bins, SCHEMA_VERSION);
	}

	/** Writes a manifest at an explicit schema version and a hash sidecar from bytes read back
	 *  from disk. This method does NOT encode anything itself -- {@code schemaVersion} only
	 *  controls which header value is written; the caller is responsible for having already
	 *  {@link IdentifierCodec#encode}d every Bin's contig_id/instance_id components before
	 *  constructing them, when writing {@link #SCHEMA_VERSION_ENCODED}. */
	public static String write(String file, String sourcePath, String sourceSha256, ArrayList<Bin> bins, int schemaVersion){
		return write(file,sourcePath,sourceSha256,"-","-",bins,schemaVersion);
	}

	/** Writes a manifest with an exact secondary-cache binding. Existing callers
	 *  without a secondary source retain byte-identical output through the
	 *  five-argument overload above. */
	public static String write(String file, String sourcePath, String sourceSha256,
			String secondarySourcePath, String secondarySourceSha256,
			ArrayList<Bin> bins, int schemaVersion){
		return write(file,sourcePath,sourceSha256,secondarySourcePath,secondarySourceSha256,bins,schemaVersion,0);
	}

	/** Writes a loaded manifest while retaining its declared split policy and source bindings. */
	public static String write(String file,Manifest manifest){
		if(manifest==null){throw new IllegalArgumentException("Cannot write a null bin manifest");}
		return write(file,manifest.sourcePath,manifest.sourceSha256,manifest.secondarySourcePath,
			manifest.secondarySourceSha256,manifest.bins,manifest.schemaVersion,manifest.splitModulus);
	}

	/** Shared implementation keeps all legacy overloads byte-compatible when no policy exists. */
	private static String write(String file,String sourcePath,String sourceSha256,
			String secondarySourcePath,String secondarySourceSha256,ArrayList<Bin> bins,int schemaVersion,int splitModulus){
		final Manifest manifest=new Manifest(sourcePath,sourceSha256,secondarySourcePath,secondarySourceSha256,bins,schemaVersion,splitModulus);
		final ByteBuilder header=formatHeader(manifest);
		final HashSet<String> ids=new HashSet<String>();
		for(Bin b : bins){
			if(!ids.add(b.id)){throw new RuntimeException("Duplicate bin_id "+b.id);}
			// Row/file version mismatch (Yoimiya source-review finding, 2026-09-09): a v2 Bin
			// written under a v1 header would emit encoded bytes that a v1 reader treats as
			// literal text (silent identity change on reload); a v1 Bin written under a v2 header
			// would have its RAW bytes wrongly decode()d on reload (a literal "%41" corrupted to
			// "A"). Reject before opening any output, not after.
			if(b.schemaVersion!=schemaVersion){
				throw new RuntimeException("Bin "+b.id+" has schemaVersion "+b.schemaVersion
					+" but write() was called with schemaVersion "+schemaVersion+" -- refusing to "
					+"emit a manifest whose declared version does not match its rows' actual encoding");
			}
			b.validateSelections();
		}
		try(Writer writer=new Writer(file,schemaVersion,header)){
			for(Bin bin : bins){writer.add(bin);}
		}
		final String digest=sha256File(file);
		final FileFormat hf=FileFormat.testOutput(file+".sha256", FileFormat.TXT, null, false, true, false, false);
		final ByteStreamWriter hsw=new ByteStreamWriter(hf); hsw.start(); hsw.print(new ByteBuilder().append(digest).append("  ").append(file).nl());
		if(hsw.poisonAndWait()){throw new RuntimeException("I/O error writing "+file+".sha256");}
		return digest;
	}

	/**
	 * Bounded-memory output counterpart to Reader. Writes bins immediately, preserving
	 * the supplied source/policy headers. It neither retains a global ID set nor emits
	 * a hash sidecar: merge callers own uniqueness checks and staging publication.
	 * The existing write() overloads retain their global preflight and sidecars.
	 */
	public static final class Writer implements AutoCloseable {

		/** Validates header metadata before opening the output; header.bins is not traversed. */
		public Writer(String file,Manifest header){this(file,header.schemaVersion,formatHeader(header));}

		private Writer(String file,int schemaVersion,ByteBuilder header){
			this.file=file; this.schemaVersion=schemaVersion;
			final FileFormat format=FileFormat.testOutput(file,FileFormat.TXT,null,false,true,false,false);
			output=new ByteStreamWriter(format); output.start(); output.print(header);
		}

		/** Emits one already validated immutable bin without retaining its contig lists. */
		public void add(Bin bin){
			if(closed){throw new IllegalStateException("Manifest writer is closed: "+file);}
			if(bin==null || bin.schemaVersion!=schemaVersion){throw new IllegalArgumentException("Bin schema does not match streaming manifest: "+file);}
			append(bin,row); output.print(row); row.clear(); rowsWritten++;
		}

		/** Number of submitted rows, valid as an output count only after close succeeds. */
		public long rowsWritten(){return rowsWritten;}

		/** Drains the asynchronous writer and fails loudly on write/flush errors. */
		@Override public void close(){
			if(!closed){
				closed=true;
				if(output.poisonAndWait()){throw new RuntimeException("I/O error writing "+file);}
			}
		}

		private final String file;
		private final int schemaVersion;
		private final ByteStreamWriter output;
		private final ByteBuilder row=new ByteBuilder(512);
		private long rowsWritten;
		private boolean closed;
	}

	/** Shared header formatter keeps streaming and complete-manifest output byte-identical. */
	private static ByteBuilder formatHeader(Manifest m){
		if(m==null){throw new IllegalArgumentException("Manifest header is null");}
		if(m.schemaVersion!=SCHEMA_VERSION && m.schemaVersion!=SCHEMA_VERSION_ENCODED){
			throw new RuntimeException("Unrecognized schema version to write: "+m.schemaVersion);
		}
		if(m.splitModulus!=0){MagQCBinSplit.validate(m.splitModulus);}
		final boolean secondary=m.secondarySourcePath!=null && !m.secondarySourcePath.equals("-");
		if(secondary!=(m.secondarySourceSha256!=null && !m.secondarySourceSha256.equals("-"))){
			throw new RuntimeException("Secondary source path/hash must both be set or both be '-'.");
		}
		if(secondary && !m.secondarySourceSha256.matches("[0-9a-f]{64}")){
			throw new RuntimeException("Secondary source SHA-256 must be 64 lowercase hexadecimal characters.");
		}
		final ByteBuilder header=new ByteBuilder().append("#schema_version\t").append(m.schemaVersion).nl()
			.append("#source\t").append(m.sourcePath==null ? "-" : m.sourcePath).append("\t#\tsha256\t")
			.append(m.sourceSha256==null ? "-" : m.sourceSha256).nl();
		if(secondary){header.append("#source_secondary\t").append(m.secondarySourcePath)
			.append("\t#\tsha256\t").append(m.secondarySourceSha256).nl();}
		if(m.splitModulus!=0){header.append("#realization_split\t").append(MagQCBinSplit.VERSION)
			.append("\tmodulus=").append(m.splitModulus).append("\tvalidation_bucket=0\n");}
		return header.append(HEADER).nl();
	}

	/** The same row parser serves sequential manifests and private indexed replay spools. */
	static Bin parseBin(LineParser1 parser,int schemaVersion,String file){
		if(parser.terms()!=11){throw new RuntimeException("Manifest row has "+parser.terms()+" fields, expected 11 in "+file);}
		// Requested fractions may use scientific notation, so retain JDK floating parsing.
		return new Bin(parser.parseString(0),parser.parseString(1),parser.parseInt(2),parser.parseString(3),
			Double.parseDouble(parser.parseString(4)),Double.parseDouble(parser.parseString(5)),
			parser.parseString(6),parser.parseString(7),parser.parseString(8),parser.parseString(9),
			parser.parseString(10),schemaVersion);
	}

	/** Appends one complete row without changing selection or membership strings. */
	static void append(Bin b, ByteBuilder out){
		out.append(b.id).tab().append(b.split).tab().append(b.targetTid).tab().append(b.contaminantTids).tab();
		out.append(Double.toString(b.compFraction)).tab().append(Double.toString(b.contFraction)).tab().append(b.spikeClass).tab();
		out.append(b.breakpoints).tab().append(b.nativeContigs).tab().append(b.foreignContigs).tab().append(b.modelRows).nl();
	}
	private static String requireToken(String s,String field){
		if(s==null || s.length()==0 || s.charAt(0)=='#' || s.indexOf('\t')>=0 || s.indexOf(';')>=0
			|| s.indexOf('\n')>=0 || s.indexOf('\r')>=0){throw new RuntimeException("Invalid "+field+": "+s);}
		return s;
	}
	private static String requireSplit(String s){
		if(!"train".equals(s) && !"val".equals(s)){throw new RuntimeException("split must be train or val, got "+s);}
		return s;
	}
	private static double checkFraction(double v,String field){
		if(Double.isNaN(v) || Double.isInfinite(v) || v<0 || v>1){throw new RuntimeException(field+" outside [0,1]: "+v);}
		return v;
	}
	private static String normalizeList(String s,String field,boolean allowDash){
		if(s==null || s.length()==0 || (allowDash && s.equals("-"))){return "-";}
		if(s.indexOf('\t')>=0 || s.indexOf(',')>=0 || s.indexOf('\n')>=0 || s.indexOf('\r')>=0){throw new RuntimeException("Invalid delimiter in "+field+": "+s);}
		for(String token : splitList(s)){
			if(token.length()==0 || token.indexOf('\t')>=0){throw new RuntimeException("Empty/malformed token in "+field);}
		}
		return s;
	}
	private static String[] splitList(String s){return "-".equals(s) ? new String[0] : s.split(";",-1);}
	private static int parseTid(String s){
		try{return Integer.parseInt(s);}catch(NumberFormatException e){throw new RuntimeException("Invalid tid: "+s);}
	}
	public static String sha256File(String file){
		try{
			final MessageDigest md=MessageDigest.getInstance("SHA-256");
			final byte[] buf=new byte[1<<16];
			try(BufferedInputStream in=new BufferedInputStream(new FileInputStream(file))){
				for(int n; (n=in.read(buf))>=0;){if(n>0){md.update(buf,0,n);}}
			}
			final StringBuilder sb=new StringBuilder(64);
			for(byte b : md.digest()){sb.append(String.format("%02x", b));}
			return sb.toString();
		}catch(IOException|NoSuchAlgorithmException e){throw new RuntimeException("Could not hash "+file,e);}
	}

	/** Small contract self-test used by CI and by the implementation handoff. */
	public static void selftest(){
		final String dir=System.getProperty("java.io.tmpdir")+"/magqc_bin_manifest_"+System.nanoTime();
		new java.io.File(dir).mkdirs();
		final ArrayList<Bin> rows=new ArrayList<Bin>();
		rows.add(new Bin("bin000000000001","train",10,"20",0.75,0.10,"ordinary","-","10|c1|n000000;10|c2|n000001","20|d1|f000000","subnetA:0;composite:4"));
		rows.add(new Bin("bin000000000002","val",20,"10",1.0,0.0,"breakpoint","g7","20|d2|n000000","-","subnetA:1"));
		final String file=dir+"/bins.tsv";
		final String h=write(file,"cache.tsv","cachehash",rows);
		final Manifest m=load(file);
		if(m.size()!=2 || !"cachehash".equals(m.sourceSha256) || !m.bins.get(1).split.equals("val")
			|| decodeSelections(m.bins.get(0).foreignContigs).size()!=1){throw new RuntimeException("Manifest round-trip failed");}
		boolean rejected=false;
		try{new Bin("bin000000000003","train",10,"20",0.5,0.1,"ordinary","-","10|c1|n000000","10|d1|f000000","x");}
		catch(RuntimeException expected){ rejected=true; }
		if(!rejected){throw new RuntimeException("Expected foreign role validation failure");}
		System.err.println("MagQCBinManifest selftest: PASS (2 rows, explicit native/foreign roles, split replay, hash)");
	}
}
