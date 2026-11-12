package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.ConfigInput;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.LongList;

/** Binds saved combined CallGenes output to its exact original panel metadata.
 * Derived from A_SampleByteFile. Streams the original once; retains every generic
 * RNA row, including unknown models, excluded families and duplicates. Other
 * feature types remain in the immutable original and explicit inventory counts.
 * No model attribution, calling, sequence reads or biological reference decisions.
 *
 * Per-input RNA builders and primitive source-line lists permit interleaved input.
 * Their aggregate backing capacities, including spare capacity, are bounded before
 * growth. Output transfers each completed builder to the native writer without a
 * payload copy. The JVM heap separately bounds maps and native I/O infrastructure.
 * Original/canonical FASTA bindings are metadata declarations, not sequence proof.
 * @author Brian Bushnell, Raiden
 */
public final class BindNcrnaCalls{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, BindNcrnaCalls.class, false);
		try{new BindNcrnaCalls(pp.args, pp.outstream).process(timer);}
		finally{Shared.closeStream(pp.outstream);}
	}
	private BindNcrnaCalls(final String[] args, final PrintStream stream){
		outstream=stream;
		Shared.setThreads(1);
		final Parser common=new Parser();
		common.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("genomes")){genomesPath=b;}
			else if(a.equals("panelcontigs")){panelContigs=b;}
			else if(a.equals("observed")){observed=b;}
			else if(a.equals("census")){census=b;}
			else if(a.equals("memberships")){memberships=b;}
			else if(a.equals("callerconfig")){callerConfig=b;}
			else if(a.equals("baseline")){baseline=b;}
			else if(a.equals("inputsha80")){inputSha80=b;}
			else if(a.equals("maxbuffer")){maxBuffer=Parse.parseKMG(b);}
			else if(a.equals("maxlines")){maxLines=Parse.parseKMG(b);}
			else if(a.equals("maxoutput")){maxOutput=Parse.parseKMG(b);}
			else if(a.equals("maxfiles")){maxFiles=Parse.parseIntKMG(b);}
			else if(!common.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		input=common.in1;
		require(input!=null && genomesPath!=null && panelContigs!=null && observed!=null && census!=null
			&& memberships!=null && callerConfig!=null && baseline!=null && common.out1!=null,
			"Supply in=combined_GFF genomes=panel_roster panelcontigs=panel_contigs observed=saved_contigs census= memberships= callerconfig= baseline= inputsha80= out=new_directory");
		require(!baseline.isEmpty() && SeedOffsetTable.validSha80(inputSha80), "Baseline and declared original-GFF sha80 are required");
		require(maxBuffer>0 && maxLines>0 && maxOutput>0 && maxFiles>METADATA_FILES, "Buffer, row, byte and file budgets must be positive");
		require(!common.overwrite && !common.append, "Call binding requires a fresh output directory");
		require(Tools.testInputFiles(false, true, input, genomesPath, panelContigs, observed, census, memberships, callerConfig), "Call-binding metadata and original GFF must be readable");
		output=new File(common.out1).getAbsoluteFile().toPath();
		require(rawSafe(output.toString()) && !Files.exists(output), "Output must be new and representable in a literal three-column query manifest");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Exact bindings       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadBindings(){
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, memberships);
		final ObjectMap<String,Genome> byPath=new ObjectMap<String,Genome>(16, String.class, Genome.class);
		for(int i=0; i<catalog.size(); i++){
			final Genome genome=catalog.get(i);
			addPath(byPath, genome.fasta, genome);
			for(int j=0; j<genome.memberships(); j++){addPath(byPath, genome.membership(j).suppliedFasta, genome);}
		}
		try(Input in=new Input(genomesPath)){
			RrnaResourceIO.header(in.reader, NcrnaCalibrationPrep.ROSTER_HEADER, genomesPath);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				parser.set(row);
				require(parser.terms()==11, "Panel genome roster requires eleven literal fields");
				final String name=text(1), assembly=text(2), fasta=text(7);
				final Genome genome=byPath.get(fasta);
				require(!name.isEmpty() && !assembly.isEmpty() && genome!=null, "Every original panel FASTA needs an exact catalog path/membership binding");
				require(parser.parseLong(9)==genome.bases && parser.parseLong(10)==genome.records, "Panel and current census dimensions differ");
				require(!byPanel.contains(name) && !byId.contains(genome.id), "Repeated panel genome or catalog input within a baseline");
				require(rawSafe(genome.id) && rawSafe(genome.fasta), "Query manifest input fields must contain no tab/CR/LF");
				final Binding binding=new Binding(name, assembly, fasta, genome, line, order.size()+1);
				byPanel.put(name, binding);
				byId.put(genome.id, binding);
				order.add(binding);
			}
		}
		require(!order.isEmpty() && (long)order.size()+METADATA_FILES<=maxFiles, "Output file budget must cover every input and diagnostic metadata file");
		try(ConfigInput config=new ConfigInput(callerConfig)){
			int index=0;
			while(config.next()){
				if(config.key.equalsIgnoreCase("in") || config.key.equalsIgnoreCase("in1")){
					require(index<order.size() && order.get(index).originalFasta.equals(config.value), "Caller config must list the exact original panel FASTAs in order");
					index++;
				}else{
					require(!config.key.equalsIgnoreCase("config") && !config.key.equalsIgnoreCase("in2"), "Nested/paired caller inputs are not covered by this flat-panel binding");
				}
			}
			require(index==order.size(), "Caller config must include the complete panel, including zero-RNA inputs");
		}
	}
	private static void addPath(final ObjectMap<String,Genome> map, final String path, final Genome genome){
		final Genome old=map.get(path);
		require(old==null || old==genome, "Declared FASTA path belongs to conflicting catalog inputs");
		map.put(path, genome);
	}

	private void loadContigs(){
		try(Input in=new Input(panelContigs)){
			final byte[] header=in.reader.nextLine();
			require(Arrays.equals(header, PANEL_CONTIGS.getBytes(StandardCharsets.US_ASCII))
				|| Arrays.equals(header, B12_CONTIGS.getBytes(StandardCharsets.US_ASCII)), "Unknown original panel-contig schema");
			Binding binding=null;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==4, "Panel contigs require genome/assembly/name/length");
				if(binding==null || !parser.termEquals(binding.name, 0)){binding=byPanel.get(text(0));}
				require(binding!=null && parser.termEquals(binding.assembly, 1), "Panel contig ownership/assembly does not match its original roster");
				probe.set(parser, 2);
				final long length=parser.parseLong(3);
				require(probe.to>probe.from && length>0 && length<=Integer.MAX_VALUE, "Panel contigs require literal names and positive int-sized lengths");
				require(!contigs.contains(probe), "A repeated literal contig name cannot attribute a combined GFF uniquely");
				final Contig contig=new Contig(probe.copy(), binding, (int)length);
				contigs.put(contig.name, contig);
				binding.contigs.add(contig);
				binding.bases=Math.addExact(binding.bases, length);
			}
		}
		for(final Binding item:order){
			require(item.contigs.size()==item.genome.records && item.bases==item.genome.bases, "Original panel contigs must conserve the bound input dimensions");
		}
		String lastId=null;
		Binding binding=null;
		try(Input in=new Input(observed)){
			RrnaResourceIO.header(in.reader, WriteNcrnaReference.CONTIG_HEADER, observed);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parser.set(row);
				require(parser.terms()==4, "Observed contigs require four fields");
				if(lastId==null || !parser.termEquals(lastId, 0)){
					lastId=text(0);
					binding=byId.get(lastId);
				}
				if(binding==null){continue;}
				require(same(1, binding.fastaEncoded) || decoded(1).equals(binding.genome.fasta),
					"Observed input FASTA differs from the exact catalog binding");
				setDecodedProbe(2);
				final Contig contig=contigs.get(probe);
				require(contig!=null && contig.owner==binding && !contig.observed && parser.parseLong(3)==contig.length,
					"Saved observed contig differs from its original panel owner/name/length");
				contig.observed=true;
			}
		}
		for(final Binding item:order){
			for(final Contig contig:item.contigs){require(contig.observed, "An original panel contig is missing from the saved observed input");}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Original call stream    ----------------*/
	/*--------------------------------------------------------------*/

	private void readCalls(){
		Contig contig=null;
		try(Input in=new Input(input)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				physicalLines++;
				if(row.length==0){
					blanks++;
					continue;
				}
				if(row[0]=='#'){
					comments++;
					continue;
				}
				parser.set(row);
				if(parser.terms()!=9 || !parser.termEquals("BBTools", 1)){
					throw new IllegalArgumentException("Expected an original nine-field CallGenes feature at physical line "+physicalLines);
				}
				final int type=parser.termEquals("CDS", 2) ? 0 : parser.termEquals("tRNA", 2) ? 1
					: parser.termEquals("rRNA", 2) ? 2 : parser.termEquals("RNA", 2) ? 3 : -1;
				require(type>=0, "Unknown feature type cannot be dropped from the original population");
				probe.set(parser, 0);
				if(contig==null || !contig.name.equals(probe)){contig=contigs.get(probe);}
				final long start=parser.parseLong(3), stop=parser.parseLong(4);
				if(contig==null || start<1 || stop<start || stop>contig.length
					|| !(parser.termEquals("+", 6) || parser.termEquals("-", 6))){
					throw new IllegalArgumentException("Original feature has unknown ownership or invalid contig bounds/strand at line "+physicalLines);
				}
				features++;
				contig.owner.features++;
				contig.owner.counts[type]++;
				totals[type]++;
				if(type==3){retain(contig.owner, row, physicalLines);}
			}
		}
		long classified=0;
		for(final long count:totals){classified+=count;}
		require(classified==features, "Every original feature must contribute once to the type inventory");
		require(physicalLines==features+comments+blanks, "Every original physical line must be a classified feature, comment or blank");
	}

	/** Capacity reservation follows ByteBuilder.expand and LongList.add growth.
	 * Both use MAX_ARRAY_LEN; check before allocating, then verify actual capacities.
	 */
	private void retain(final Binding binding, final byte[] row, final long sourceLine){
		require(totals[3]<=maxLines, "Retained RNA line budget exceeded");
		final int oldBytes=binding.rna==null ? 0 : binding.rna.array.length;
		final int oldLongs=binding.lines==null ? 0 : binding.lines.array.length;
		final int size=binding.lines==null ? 0 : binding.lines.size;
		final long needed=(binding.rna==null ? GFF_HEADER.length() : binding.rna.length())+(long)row.length+1;
		require(needed<=ByteBuilder.MAXLEN, "One input's retained RNA bytes exceed native array capacity");
		long bytes=oldBytes==0 ? 128 : oldBytes;
		if(bytes<needed){while(bytes<=needed){bytes<<=1;}}
		bytes=Math.min(bytes, ByteBuilder.MAXLEN);
		final long longs=oldLongs==0 ? 16 : size>=oldLongs ? Math.min(Shared.MAX_ARRAY_LEN, size*2L+1) : oldLongs;
		final long next=retainedCapacity-oldBytes-8L*oldLongs+bytes+8L*longs;
		require(next<=maxBuffer, "Aggregate retained RNA buffer/line-array capacity budget exceeded");
		if(binding.rna==null){
			binding.rna=new ByteBuilder((int)bytes);
			binding.rna.append(GFF_HEADER);
			binding.lines=new LongList((int)longs);
		}
		binding.rna.ensureExtra(row.length+1);
		binding.rna.append(row).nl();
		binding.lines.add(sourceLine);
		final long actual=retainedCapacity-oldBytes-8L*oldLongs+binding.rna.array.length+8L*binding.lines.array.length;
		require(actual==next && actual<=maxBuffer, "Native buffer growth changed; update the pre-allocation capacity reservation");
		retainedCapacity=actual;
		peakCapacity=Math.max(peakCapacity, actual);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Output           ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		loadBindings();
		loadContigs();
		readCalls();
		try{Files.createDirectory(output);}
		catch(IOException failure){throw new IllegalStateException("Could not create a fresh query-projection directory", failure);}
		writeSources();
		writeInputs();
		writeSummary();
		writeManifest();
		timer.stop();
		outstream.println("Bound "+totals[3]+" original RNA rows across "+order.size()+" inputs; logical output bytes="+outputBytes+". Time: "+timer);
		outstream.println("Projection only: original non-RNA calls remain in the source and original_counts.tsv; no calls or references were inferred.");
	}
	private void writeSources(){
		try(Output out=new Output(path("sources.tsv"))){
			header(out, "role\tpath\tdeclared_sha80\tbaseline\tqualification");
			final String[] roles={"original_gff", "panel_genomes", "panel_contigs", "observed_contigs", "census", "memberships", "caller_config"};
			final String[] paths={input, genomesPath, panelContigs, observed, census, memberships, callerConfig};
			for(int i=0; i<roles.length; i++){
				buffer.clear().append(roles[i]).tab();
				field(new File(paths[i]).getAbsolutePath()).tab().append(i==0 ? inputSha80 : "NA").tab();
				field(baseline).tab().append("DECLARED_METADATA_IDENTITY_NOT_SEQUENCE_VERIFICATION").nl();
				write(out, buffer);
			}
		}
	}
	private void writeInputs(){
		long copied=0;
		try(Output bindings=new Output(path("bindings.tsv")); Output counts=new Output(path("original_counts.tsv"));
			Output lines=new Output(path("line_map.tsv.gz")); Output observedOut=new Output(path("contigs.tsv.gz"))){
			header(bindings, "input_id\tpanel_genome\tassembly\toriginal_fasta\tcatalog_fasta\tpanel_line\tquery_gff\tprojection");
			header(counts, "input_id\tsource_features\tCDS\ttRNA\trRNA\tprojected_RNA\tquery_features\tqualification");
			header(lines, "source_physical_line\tinput_id\tquery_physical_line");
			header(observedOut, WriteNcrnaReference.CONTIG_HEADER);
			for(final Binding binding:order){
				final String gff=queryPath(binding);
				final long capacity=binding.rna==null ? 0 : binding.rna.array.length+8L*binding.lines.array.length;
				final ByteBuilder payload;
				if(binding.rna==null){
					require(retainedCapacity+32<=maxBuffer, "Empty query header exceeds the retained-buffer budget");
					payload=new ByteBuilder(32).append(GFF_HEADER);
					peakCapacity=Math.max(peakCapacity, retainedCapacity+32);
					emptyInputs++;
				}else{payload=binding.rna;}
				charge(payload.length());
				binding.rna=null;
				try(Output out=new Output(gff)){
					// No previous print: the native writer buffer is empty, so this
					// transfers ownership (ByteStreamWriter.print:924), without copying.
					out.writer.print(payload, true);
				}
				if(binding.lines!=null){
					require(binding.lines.size==binding.counts[3], "Every retained RNA row needs one original physical-line witness");
					for(int i=0; i<binding.lines.size; i++){
						assert(i==0 || binding.lines.get(i)>binding.lines.get(i-1)) : "Per-input projection must preserve original call order, including duplicates";
						buffer.clear().append(binding.lines.get(i)).tab();
						field(binding.genome.id).tab().append(i+2L).nl();
						write(lines, buffer);
						copied++;
					}
				}
				binding.lines=null;
				retainedCapacity-=capacity;
				buffer.clear();
				field(binding.genome.id).tab();
				field(binding.name).tab();
				field(binding.assembly).tab();
				field(binding.originalFasta).tab();
				field(binding.genome.fasta).tab().append(binding.line).tab();
				field(gff).append("\tGENERIC_RNA_ONLY_ORIGINAL_OTHER_FEATURES_RETAINED\n");
				write(bindings, buffer);
				long classified=0;
				for(final long count:binding.counts){classified+=count;}
				require(classified==binding.features, "Per-input type counts must conserve all original features");
				buffer.clear();
				field(binding.genome.id).tab().append(binding.features);
				for(final long count:binding.counts){buffer.tab().append(count);}
				buffer.tab().append(binding.counts[3]).append("\tPROJECTED_ZERO_OTHER_TYPES_IS_NOT_SOURCE_ABSENCE\n");
				write(counts, buffer);
				for(final Contig contig:binding.contigs){
					buffer.clear();
					field(binding.genome.id).tab();
					field(binding.genome.fasta).tab();
					BuildNcrnaReference.field(buffer, contig.name.bytes, 0, contig.name.bytes.length).tab().append(contig.length).nl();
					write(observedOut, buffer);
				}
			}
		}
		require(copied==totals[3] && retainedCapacity==0, "Original RNA rows, projected rows and physical-line witnesses must be bijective");
	}
	private void writeSummary(){
		try(Output out=new Output(path("summary.tsv"))){
			header(out, "baseline\tinputs\tphysical_lines\tCDS\ttRNA\trRNA\tprojected_RNA\tcomments\tblanks\tempty_RNA_inputs\tpeak_retained_array_bytes\tmax_retained_array_bytes\tmax_RNA_rows\toutput_files\tmax_output_files\tmax_logical_output_bytes\tqualification");
			buffer.clear();
			field(baseline).tab().append(order.size()).tab().append(physicalLines);
			for(final long count:totals){buffer.tab().append(count);}
			buffer.tab().append(comments).tab().append(blanks).tab().append(emptyInputs).tab().append(peakCapacity)
				.tab().append(maxBuffer).tab().append(maxLines).tab().append(order.size()+METADATA_FILES).tab().append(maxFiles)
				.tab().append(maxOutput).append("\tQUERY_PROJECTION_NOT_BIOLOGICAL_REFERENCE\n");
			write(out, buffer);
		}
	}
	/** Final query manifest follows every reader and payload close; it is never truth. */
	private void writeManifest(){
		final String pending=path("calls.partial.tsv");
		try(Output out=new Output(pending)){
			header(out, GradeNcrnaReference.MANIFEST_HEADER);
			for(final Binding binding:order){
				buffer.clear().append(binding.genome.id.getBytes(StandardCharsets.UTF_8)).tab()
					.append(binding.genome.fasta.getBytes(StandardCharsets.UTF_8)).tab().append(queryPath(binding).getBytes(StandardCharsets.UTF_8)).nl();
				write(out, buffer);
			}
		}
		try{Files.move(new File(pending).toPath(), output.resolve("calls.tsv"), StandardCopyOption.ATOMIC_MOVE);}
		catch(IOException failure){throw new IllegalStateException("Could not publish completed query projection", failure);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private void charge(final int bytes){
		require(bytes>=0 && bytes<=maxOutput-outputBytes, "Logical output-byte budget exceeded; partial directory is not complete");
		outputBytes+=bytes;
	}
	private void write(final Output out, final ByteBuilder bytes){
		charge(bytes.length());
		out.writer.print(bytes);
	}
	private void header(final Output out, final String header){write(out, buffer.clear().append(header).nl());}
	private ByteBuilder field(final String value){return BuildNcrnaReference.field(buffer, value);}
	private String path(final String name){return output.resolve(name).toString();}
	private String queryPath(final Binding binding){return path("input_"+binding.ordinal+".gff.gz");}
	private String text(final int column){
		parser.setBounds(column);
		final String result=new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		require(result.indexOf('\uFFFD')<0, "Malformed UTF-8 changes metadata identity");
		return result;
	}
	private String decoded(final int column){
		parser.setBounds(column);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	/** Compares canonical escaped bytes without constructing a path String per contig. */
	private boolean same(final int column, final byte[] expected){
		parser.setBounds(column);
		if(parser.b()-parser.a()!=expected.length){return false;}
		final byte[] line=parser.line();
		for(int i=0; i<expected.length; i++){
			if(line[parser.a()+i]!=expected[i]){return false;}
		}
		return true;
	}
	private void setDecodedProbe(final int column){
		parser.setBounds(column);
		for(int i=parser.a(); i<parser.b(); i++){
			if(parser.line()[i]=='%'){
				final byte[] bytes=decoded(column).getBytes(StandardCharsets.UTF_8);
				probe.set(bytes, 0, bytes.length);
				return;
			}
		}
		probe.set(parser, column);
	}
	private static boolean rawSafe(final String value){return value.indexOf('\t')<0 && value.indexOf('\r')<0 && value.indexOf('\n')<0;}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	private static final class Binding{
		Binding(final String name_, final String assembly_, final String fasta_, final Genome genome_, final long line_, final int ordinal_){
			name=name_;
			assembly=assembly_;
			originalFasta=fasta_;
			genome=genome_;
			fastaEncoded=BuildNcrnaReference.field(new ByteBuilder(), genome.fasta).toBytes();
			line=line_;
			ordinal=ordinal_;
		}
		final String name, assembly, originalFasta;
		final Genome genome;
		final byte[] fastaEncoded;
		final long line;
		final int ordinal;
		final ArrayList<Contig> contigs=new ArrayList<Contig>();
		final long[] counts=new long[4];
		long bases, features;
		ByteBuilder rna;
		LongList lines;
	}
	private static final class Contig{
		Contig(final Target name_, final Binding owner_, final int length_){
			name=name_;
			owner=owner_;
			length=length_;
		}
		final Target name;
		final Binding owner;
		final int length;
		boolean observed;
	}
	private String input, genomesPath, panelContigs, observed, census, memberships, callerConfig, baseline, inputSha80;
	private long maxBuffer=128L*1024*1024, maxLines=2000000, maxOutput=256L*1024*1024;
	private int maxFiles=10000;
	private final Path output;
	private final PrintStream outstream;
	private final LineParser1 parser=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final Target probe=new Target();
	private final ObjectMap<String,Binding> byPanel=new ObjectMap<String,Binding>(16, String.class, Binding.class);
	private final ObjectMap<String,Binding> byId=new ObjectMap<String,Binding>(16, String.class, Binding.class);
	private final ObjectMap<Target,Contig> contigs=new ObjectMap<Target,Contig>(16, Target.class, Contig.class);
	private final ArrayList<Binding> order=new ArrayList<Binding>();
	private final ByteBuilder buffer=new ByteBuilder();
	private final long[] totals=new long[4];
	private long physicalLines, features, comments, blanks, retainedCapacity, peakCapacity, outputBytes, emptyInputs;
	private static final int METADATA_FILES=7;
	private static final String GFF_HEADER="##gff-version 3\n";
	public static final String PANEL_CONTIGS="genome\tassembly\tcontig\tlength";
	public static final String B12_CONTIGS="genome\tassembly\tseqid\tlength";
}
