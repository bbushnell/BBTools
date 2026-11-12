package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Locale;

import map.IntHashMap;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog.Genome;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import structures.IntList;

/** Partitions whole inputs so literal contig names identify their owner within a batch.
 * Adapted from template.A_SampleByteFile. One thread reads saved metadata and owns
 * retained byte keys and primitive lists; no FASTA, calls or annotations are read.
 * Original census order defines deterministic first-fit placement. The resulting
 * batch count is an upper bound, not a proof of the minimum coloring. Maximum
 * input multiplicity of one name supplies a separately reported lower bound.
 * No names, paths, sequences, taxonomy or historical memberships are rewritten.
 * @author Brian Bushnell, Raiden
 */
public final class PartitionNcrnaBenchmarkInputs{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PartitionNcrnaBenchmarkInputs tool=new PartitionNcrnaBenchmarkInputs(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	private PartitionNcrnaBenchmarkInputs(final String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		Shared.setThreads(1);
		final Parser common=new Parser();
		common.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("census")){census=b;}
			else if(a.equals("memberships")){memberships=b;}
			else if(a.equals("contigs")){contigs=b;}
			else if(!common.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		selection=common.in1;
		require(selection!=null && census!=null && memberships!=null && contigs!=null && common.out1!=null,
			"Supply in=selection.tsv census= memberships= contigs= out=fresh_directory");
		require(!common.overwrite && !common.append && common.in2==null && common.out2==null
			&& common.maxReads<0 && common.samplerate==1, "Partition complete saved populations into fresh outputs; no read limits or paired inputs");
		require(Tools.testInputFiles(false, true, selection, census, memberships, contigs), "All partition metadata must be readable");
		root=new File(common.out1).getAbsoluteFile();
		require(!root.exists() && root.getParentFile().isDirectory(), "Partition output needs a fresh directory under an existing parent");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Saved metadata        ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		loadSelection();
		readContigs();
		partition();
		try{Files.createDirectory(root.toPath());}
		catch(IOException failure){throw new IllegalArgumentException("Cannot create partition output "+root, failure);}
		write();
		timer.stop();
		outstream.println("Partitioned "+selected.size+" inputs into "+batches.size()+" noncolliding batches; multiplicity lower bound="+lowerBound+". Time: "+timer);
		outstream.println("NCRNA_BENCHMARK_PARTITION_COMPLETE; saved metadata only; greedy count is not a proven minimum.");
	}

	private void loadSelection(){
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, memberships);
		inputs=new Member[catalog.size()];
		for(int i=0; i<inputs.length; i++){
			final Member member=new Member(catalog.get(i), i);
			inputs[i]=member;
			byId.put(member.genome.id, member);
		}
		int last=-1;
		try(Input in=new Input(selection)){
			RrnaResourceIO.header(in.reader, SELECTION_HEADER, selection);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				parser.set(row);
				require(parser.terms()==8, "Use the complete eight-column native selection table");
				final Member member=byId.get(parser.parseString(0));
				require(member!=null && member.index>last, "Selection must contain unique catalog IDs in original census order");
				require(same(1, member.fasta) && parser.parseLong(2)==member.genome.bases
					&& parser.parseLong(3)==member.genome.records, "Selected path and dimensions must match the exact catalog input");
				member.selected=true;
				selected.add(member.index);
				last=member.index;
			}
		}
		require(selected.size>0, "An empty selection cannot define grading batches");
	}

	/** Complete catalog dimensions are checked even when only a subset is selected.
	 * Worker-buffered observations may interleave inputs; no block order is assumed.
	 * Only selected names are retained and checked for within-input duplication.
	 */
	private void readContigs(){
		Member current=null;
		try(Input in=new Input(contigs)){
			RrnaResourceIO.header(in.reader, WriteNcrnaReference.CONTIG_HEADER, contigs);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				parser.set(row);
				require(parser.terms()==4, "Saved contigs require input ID, escaped FASTA, escaped literal name and length");
				if(current==null || !parser.termEquals(current.genome.id, 0)){
					current=byId.get(parser.parseString(0));
					require(current!=null, "Saved contigs must identify an original catalog input");
					current.seen=true;
				}
				require(same(1, current.encodedFasta) || decoded(1).equals(current.genome.fasta),
					"Saved contig FASTA must retain its exact catalog binding");
				final long length=parser.parseLong(3);
				require(length>0 && length<=Integer.MAX_VALUE, "Contig length must be positive and representable by the native binder");
				if(current.selected){
					decodeName(2);
					final int name=internName();
					if(current.seenNames==null){current.seenNames=new IntHashMap(1);}
					require(current.seenNames.get(name)<0, "A literal contig name repeats within input "+current.genome.id);
					current.seenNames.put(name, 1);
					multiplicity.set(name, Math.addExact(multiplicity.get(name), 1));
					current.names.add(name);
					current.lengths.add((int)length);
				}
				current.records=Math.addExact(current.records, 1);
				current.bases=Math.addExact(current.bases, length);
				require(current.records<=current.genome.records && current.bases<=current.genome.bases,
					"Saved contig dimensions exceed catalog input "+current.genome.id);
			}
		}
		for(final Member member:inputs){
			require(member.seen, "Saved contigs are missing catalog input "+member.genome.id);
			checkDimensions(member);
			assert(!member.selected || member.names.size==member.records && member.lengths.size==member.names.size)
				: "Batch writers need one aligned name/length entry for every selected input observation";
			member.seenNames=null;
		}
	}

	private static void checkDimensions(final Member member){
		require(member.records==member.genome.records && member.bases==member.genome.bases,
			"Saved contig dimensions disagree with catalog input "+member.genome.id);
	}

	/** Unescaped bytes are borrowed only until lookup; retained dictionary keys are copied.
	 * Escaped identifiers use the same strict UTF-8 decoder as the existing binder.
	 */
	private void decodeName(final int field){
		parser.setBounds(field);
		probe.set(parser, field);
		for(int i=parser.a(); i<parser.b(); i++){
			if(parser.line()[i]=='%'){
				final byte[] decoded=decoder.decode(parser.line(), parser.a(), parser.b()).getBytes(StandardCharsets.UTF_8);
				probe.set(decoded, 0, decoded.length);
				break;
			}
		}
		require(probe.to>probe.from, "A literal contig identifier cannot be empty");
		for(int i=probe.from; i<probe.to; i++){
			require(probe.bytes[i]!='\t' && probe.bytes[i]!='\r' && probe.bytes[i]!='\n',
				"Decoded contig names cannot contain tab, CR or LF in a literal binder roster");
		}
	}

	private int internName(){
		final int old=nameIds.get(probe);
		if(old>=0){return old;}
		final int index=names.size();
		final Target retained=probe.copy();
		nameIds.put(retained, index);
		names.add(retained);
		multiplicity.add(0);
		assert(nameIds.size()==names.size()) : "Each retained contig key owns one primitive name index";
		return index;
	}

	/*--------------------------------------------------------------*/
	/*----------------       First-fit batches      ----------------*/
	/*--------------------------------------------------------------*/

	private void partition(){
		final IntList[] used=new IntList[names.size()];
		for(int n=0; n<names.size(); n++){
			final int count=multiplicity.get(n);
			lowerBound=Math.max(lowerBound, count);
			if(count>1){
				sharedNames++;
				sharedOccurrences=Math.addExact(sharedOccurrences, count);
				used[n]=new IntList(Math.min(count, 4));
			}
		}
		final int[] forbidden=new int[selected.size];
		for(int ordinal=0; ordinal<selected.size; ordinal++){
			final Member member=inputs[selected.get(ordinal)];
			final int epoch=ordinal+1;
			boolean shared=false;
			for(int i=0; i<member.names.size; i++){
				final IntList previous=used[member.names.get(i)];
				if(previous!=null){
					shared=true;
					for(int j=0; j<previous.size; j++){forbidden[previous.get(j)]=epoch;}
				}
			}
			if(shared){sharedInputs++;}
			int batch=0;
			while(batch<batches.size() && forbidden[batch]==epoch){batch++;}
			if(batch==batches.size()){batches.add(new IntList(4));}
			member.batch=batch;
			member.local=batches.get(batch).size+1;
			batches.get(batch).add(member.index);
			for(int i=0; i<member.names.size; i++){
				final IntList previous=used[member.names.get(i)];
				if(previous!=null){previous.add(batch);}
			}
		}
		assert(batches.size()>=lowerBound && batches.size()<=selected.size)
			: "Every shared name needs distinct batches, and assigning one input per batch is always feasible";
	}

	/*--------------------------------------------------------------*/
	/*----------------       Binder-ready outputs   ----------------*/
	/*--------------------------------------------------------------*/

	private void write(){
		try(Output assignments=new Output(path("assignments.tsv"))){
			assignments.writer.println("input_id\tfasta\tbatch\tbatch_input_index");
			for(int i=0; i<selected.size; i++){
				final Member member=inputs[selected.get(i)];
				buffer.clear().append(member.genome.id).tab().append(member.fasta).tab().append(member.batch).tab().append(member.local).nl();
				assignments.writer.print(buffer);
			}
		}
		try(Output batchTotals=new Output(path("batches.tsv"))){
			batchTotals.writer.println("batch\tinputs\tnucleotide_bases\tfasta_records\tbinder_maxfiles");
			for(int b=0; b<batches.size(); b++){writeBatch(b, batchTotals);}
		}
		long bases=0, records=0;
		for(int i=0; i<selected.size; i++){
			final Genome genome=inputs[selected.get(i)].genome;
			bases=Math.addExact(bases, genome.bases);
			records=Math.addExact(records, genome.records);
		}
		try(Output summary=new Output(path("summary.tsv"))){
			summary.writer.println("selected_inputs\tnucleotide_bases\tfasta_records\tdistinct_contig_names\tshared_contig_names\tinputs_with_shared_names\tshared_name_occurrences\tmax_input_multiplicity\tgreedy_batches\tqualification");
			buffer.clear().append(selected.size).tab().append(bases).tab().append(records).tab().append(names.size())
				.tab().append(sharedNames).tab().append(sharedInputs).tab().append(sharedOccurrences).tab().append(lowerBound).tab().append(batches.size())
				.append("\tSaved metadata; greedy batch count is not a proven minimum\n");
			summary.writer.print(buffer);
		}
	}

	private void writeBatch(final int batch, final Output batchTotals){
		final IntList list=batches.get(batch);
		long bases=0, records=0;
		try(Output config=new Output(path("batch_"+batch+".inputs.config")); Output genomes=new Output(path("batch_"+batch+".genomes.tsv"));
			Output targets=new Output(path("batch_"+batch+".contigs.tsv")); Output limits=new Output(path("batch_"+batch+".binder.config"))){
			genomes.writer.println(NcrnaCalibrationPrep.ROSTER_HEADER);
			targets.writer.println(BindNcrnaCalls.PANEL_CONTIGS);
			for(int i=0; i<list.size; i++){
				final Member member=inputs[list.get(i)];
				bases=Math.addExact(bases, member.genome.bases);
				records=Math.addExact(records, member.genome.records);
				buffer.clear().append("in=").append(member.fasta).nl();
				config.writer.print(buffer);
				buffer.clear().append(i+1).tab().append(member.genome.id).append("\tNA\tNA\tNA\tNA\tNA\t")
					.append(member.fasta).append("\tNA\t").append(member.genome.bases).tab().append(member.genome.records).nl();
				genomes.writer.print(buffer);
				for(int n=0; n<member.names.size; n++){
					buffer.clear().append(member.genome.id).append("\tNA\t").append(names.get(member.names.get(n)).bytes)
						.tab().append(member.lengths.get(n)).nl();
					targets.writer.print(buffer);
				}
			}
			buffer.clear().append("maxfiles=").append(Math.addExact(list.size, 7)).nl();
			limits.writer.print(buffer);
		}
		buffer.clear().append(batch).tab().append(list.size).tab().append(bases).tab().append(records)
			.tab().append(Math.addExact(list.size, 7)).nl();
		batchTotals.writer.print(buffer);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Helpers and state      ----------------*/
	/*--------------------------------------------------------------*/

	private boolean same(final int field, final byte[] value){
		parser.setBounds(field);
		if(parser.b()-parser.a()!=value.length){return false;}
		for(int i=0; i<value.length; i++){if(parser.line()[parser.a()+i]!=value[i]){return false;}}
		return true;
	}
	private String decoded(final int field){
		parser.setBounds(field);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	private String path(final String name){return new File(root, name).toString();}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	private static final class Member{
		Member(final Genome genome_, final int index_){
			genome=genome_;
			index=index_;
			require(genome.fasta.equals(genome.fasta.trim()), "FASTA paths cannot have whitespace trimmed by native config readers");
			fasta=genome.fasta.getBytes(StandardCharsets.UTF_8);
			for(final byte b:fasta){require(b!='\t' && b!='\r' && b!='\n' && b!=',', "FASTA paths must fit literal caller configs and binder rosters");}
			encodedFasta=BuildNcrnaReference.field(new ByteBuilder(), genome.fasta).toBytes();
		}
		final Genome genome;
		final int index;
		final byte[] fasta, encodedFasta;
		final IntList names=new IntList(1), lengths=new IntList(1);
		IntHashMap seenNames;
		boolean selected, seen;
		long bases, records;
		int batch, local;
	}

	private String selection, census, memberships, contigs;
	private final File root;
	private final PrintStream outstream;
	private Member[] inputs;
	private int lowerBound, sharedNames, sharedInputs;
	private long sharedOccurrences;
	private final ObjectMap<String,Member> byId=new ObjectMap<String,Member>(String.class, Member.class);
	private final ObjectIntMap<Target> nameIds=new ObjectIntMap<Target>(Target.class);
	private final ArrayList<Target> names=new ArrayList<Target>();
	private final ArrayList<IntList> batches=new ArrayList<IntList>();
	private final IntList selected=new IntList(), multiplicity=new IntList();
	private final LineParser1 parser=new LineParser1('\t');
	private final Target probe=new Target();
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final ByteBuilder buffer=new ByteBuilder();
	private static final String SELECTION_HEADER="input_id\tfasta\tnucleotide_bases\tfasta_records\tdomain\tphylum\tbase_log2\trecord_log2";
}
