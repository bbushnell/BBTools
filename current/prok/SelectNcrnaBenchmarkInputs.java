package prok;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;
import java.util.Random;

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

/** Family-blind benchmark sampling from a saved whole-input census and taxonomy.
 * Adapted from template.A_SampleByteFile; one thread owns native byte streams and
 * retained metadata. No sequences, annotations, reference counts or calls are read.
 * Each domain/phylum receives one input, then remaining slots are proportional
 * to its remaining population using largest remainders. Unknown labels survive.
 * Within each group, systematic sampling of joint log2(base-count)/log2(contig-count)
 * bins preserves workload variety. Seeded retained ranks break ties within bins.
 * Output order remains the original census order. A full-size selection includes
 * every original path. Sampling describes the supplied corpus, not all prokaryotes.
 * @author Brian Bushnell, Raiden
 */
public final class SelectNcrnaBenchmarkInputs{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final SelectNcrnaBenchmarkInputs tool=new SelectNcrnaBenchmarkInputs(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	private SelectNcrnaBenchmarkInputs(final String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		Shared.setThreads(1);
		final Parser common=new Parser();
		common.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("memberships")){memberships=b;}
			else if(a.equals("taxonomy")){taxonomy=b;}
			else if(a.equals("sizes")){sizes=Parse.parseIntArray(b, ",");}
			else if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(!common.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		census=common.in1;
		prefix=common.out1;
		require(census!=null && memberships!=null && taxonomy!=null && prefix!=null && sizes!=null && sizes.length>0,
			"Supply in=census memberships= taxonomy= sizes= out=fresh_prefix");
		require(!common.overwrite && !common.append, "Benchmark selections require fresh outputs");
		require(common.maxReads<0 && common.samplerate==1 && common.in2==null && common.out2==null,
			"Selection sizes, rather than sequence read limits or paired inputs, define these metadata populations");
		require(Tools.testInputFiles(false, true, census, memberships, taxonomy), "All saved metadata inputs must be readable");
		for(int i=0; i<sizes.length; i++){
			require(sizes[i]>0, "Selection sizes must be positive");
			for(int j=0; j<i; j++){require(sizes[i]!=sizes[j], "Repeated selection size");}
			RrnaResourceIO.fresh(prefix+"."+sizes[i], SUFFIXES);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved metadata         ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, memberships);
		members=new Member[catalog.size()];
		final Random random=new Random(seed);
		for(int i=0; i<members.length; i++){
			final Member member=new Member(catalog.get(i), i, random.nextLong());
			members[i]=member;
			final Target id=new Target();
			id.set(member.id, 0, member.id.length);
			ids.put(id, i);
		}
		assert(ids.size()==members.length) : "The catalog's unique input IDs must retain one dense selection index each";
		readTaxonomy();
		groups.sort(GROUP_ORDER);
		for(final Group group:groups){
			group.sorted=group.members.toArray(new Member[group.members.size()]);
			Arrays.sort(group.sorted, MEMBER_ORDER);
		}
		for(final int size:sizes){select(size);}
		timer.stop();
		outstream.println("Selected from "+members.length+" catalog inputs in "+groups.size()+" domain/phylum groups. Time: "+timer);
		outstream.println("NCRNA_BENCHMARK_INPUT_SELECTION_COMPLETE; family-blind metadata only; corpus scope remains declared.");
	}

	/** Exact IDs and FASTA bytes must cover the whole census once; taxonomy is not guessed. */
	private void readTaxonomy(){
		try(Input in=new Input(taxonomy)){
			RrnaResourceIO.header(in.reader, TAXONOMY_HEADER, taxonomy);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				parser.set(row);
				require(parser.terms()==7, "Taxonomy requires the complete seven-field native schema");
				probe.set(parser, 0);
				final int index=ids.get(probe);
				require(index>=0, "Taxonomy names an input outside the saved census");
				final Member member=members[index];
				probe.set(parser, 1);
				require(member.group==null && equal(probe, member.fasta), "Taxonomy must bind each exact input/path once");
				key.clear();
				appendField(key, 2).tab();
				appendField(key, 3);
				probe.set(key.array, 0, key.length());
				Group group=byGroup.get(probe);
				if(group==null){
					final byte[] domain=copyField(2), phylum=copyField(3);
					require(domain.length>0 && phylum.length>0, "Missing taxonomy must use explicit labels, not empty fields");
					group=new Group(domain, phylum);
					byGroup.put(probe.copy(), group);
					groups.add(group);
				}
				group.members.add(member);
				member.group=group;
			}
		}
		for(final Member member:members){require(member.group!=null, "Every census input needs explicit taxonomy, including unknowns");}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Proportional selection   ----------------*/
	/*--------------------------------------------------------------*/

	private void select(final int size){
		require(size>=groups.size() && size<=members.length,
			"Selection must retain every domain/phylum group and cannot exceed the complete census");
		final int slots=size-groups.size(), population=members.length-groups.size();
		int allocated=groups.size();
		for(final Group group:groups){
			final long product=Math.multiplyExact((long)slots, group.sorted.length-1);
			group.quota=1+(population==0 ? 0 : (int)(product/population));
			group.remainder=population==0 ? 0 : product%population;
			allocated+=group.quota-1;
		}
		final Group[] remainders=groups.toArray(new Group[groups.size()]);
		Arrays.sort(remainders, REMAINDER_ORDER);
		for(int i=0; allocated<size; i++, allocated++){
			require(i<remainders.length && remainders[i].quota<remainders[i].sorted.length,
				"Largest remainders must distribute only unallocated population slots");
			remainders[i].quota++;
		}
		final boolean[] chosen=new boolean[members.length];
		int total=0;
		for(final Group group:groups){
			for(int i=0; i<group.quota; i++){
				final int position=(int)(((2L*i+1)*group.sorted.length)/(2L*group.quota));
				final Member member=group.sorted[position];
				require(!chosen[member.index], "Systematic sampling must select each original input at most once");
				chosen[member.index]=true;
				total++;
			}
		}
		require(total==size, "Allocated strata must conserve the exact requested input count");
		write(size, chosen);
	}

	/** Writes only retained metadata and native caller input declarations. */
	private void write(final int size, final boolean[] chosen){
		final String out=prefix+"."+size;
		try(Output list=new Output(out+".inputs.txt"); Output config=new Output(out+".inputs.config");
			Output flat=new Output(out+".flat.config"); Output selection=new Output(out+".selection.tsv");
			Output composition=new Output(out+".composition.tsv"); Output totals=new Output(out+".totals.tsv")){
			selection.writer.println("input_id\tfasta\tnucleotide_bases\tfasta_records\tdomain\tphylum\tbase_log2\trecord_log2");
			config.writer.print("in=");
			boolean first=true;
			long bases=0, records=0, selectedBases=0, selectedRecords=0;
			for(final Member member:members){
				bases=Math.addExact(bases, member.genome.bases);
				records=Math.addExact(records, member.genome.records);
				if(!chosen[member.index]){continue;}
				selectedBases=Math.addExact(selectedBases, member.genome.bases);
				selectedRecords=Math.addExact(selectedRecords, member.genome.records);
				list.writer.println(member.fasta);
				buffer.clear().append("in=").append(member.fasta).nl();
				flat.writer.print(buffer);
				buffer.clear();
				if(!first){buffer.append(',');}
				first=false;
				config.writer.print(buffer.append(member.fasta));
				buffer.clear().append(member.id).tab().append(member.fasta).tab().append(member.genome.bases)
					.tab().append(member.genome.records).tab().append(member.group.domain).tab().append(member.group.phylum)
					.tab().append(member.baseBin).tab().append(member.recordBin).nl();
				selection.writer.print(buffer);
			}
			config.writer.print("\n");
			composition.writer.println("domain\tphylum\tbase_log2\trecord_log2\tcorpus_inputs\tselected_inputs\tcorpus_bases\tselected_bases\tcorpus_records\tselected_records");
			for(final Group group:groups){writeComposition(group, chosen, composition);}
			totals.writer.println("corpus_inputs\tselected_inputs\tcorpus_bases\tselected_bases\tcorpus_records\tselected_records");
			buffer.clear().append(members.length).tab().append(size).tab().append(bases).tab().append(selectedBases)
				.tab().append(records).tab().append(selectedRecords).nl();
			totals.writer.print(buffer);
		}
	}

	private void writeComposition(final Group group, final boolean[] chosen, final Output out){
		for(int from=0; from<group.sorted.length;){
			final Member first=group.sorted[from];
			int to=from, selected=0;
			long bases=0, selectedBases=0, records=0, selectedRecords=0;
			while(to<group.sorted.length && group.sorted[to].baseBin==first.baseBin && group.sorted[to].recordBin==first.recordBin){
				final Member member=group.sorted[to++];
				bases=Math.addExact(bases, member.genome.bases);
				records=Math.addExact(records, member.genome.records);
				if(chosen[member.index]){
					selected++;
					selectedBases=Math.addExact(selectedBases, member.genome.bases);
					selectedRecords=Math.addExact(selectedRecords, member.genome.records);
				}
			}
			buffer.clear().append(group.domain).tab().append(group.phylum).tab().append(first.baseBin).tab().append(first.recordBin)
				.tab().append(to-from).tab().append(selected).tab().append(bases).tab().append(selectedBases)
				.tab().append(records).tab().append(selectedRecords).nl();
			out.writer.print(buffer);
			from=to;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private ByteBuilder appendField(final ByteBuilder out, final int field){
		parser.setBounds(field);
		return out.append(parser.line(), parser.a(), parser.b()-parser.a());
	}
	private byte[] copyField(final int field){
		parser.setBounds(field);
		return Arrays.copyOfRange(parser.line(), parser.a(), parser.b());
	}
	private static boolean equal(final Target key, final byte[] bytes){
		if(key.to-key.from!=bytes.length){return false;}
		for(int i=0; i<bytes.length; i++){if(key.bytes[key.from+i]!=bytes[i]){return false;}}
		return true;
	}
	private static int compare(final byte[] a, final byte[] b){
		for(int i=0, n=Math.min(a.length, b.length); i<n; i++){
			final int c=Integer.compare(a[i]&255, b[i]&255);
			if(c!=0){return c;}
		}
		return Integer.compare(a.length, b.length);
	}
	private static void require(final boolean condition, final String message){if(!condition){throw new IllegalArgumentException(message);}}

	private static final class Member{
		Member(final Genome genome_, final int index_, final long rank_){
			genome=genome_;
			index=index_;
			rank=rank_;
			id=genome.id.getBytes(StandardCharsets.US_ASCII);
			fasta=genome.fasta.getBytes(StandardCharsets.UTF_8);
			require(genome.fasta.indexOf(',')<0 && genome.fasta.indexOf('\t')<0 && genome.fasta.indexOf('\n')<0 && genome.fasta.indexOf('\r')<0,
				"Caller input declarations require paths without comma or record delimiters");
			baseBin=63-Long.numberOfLeadingZeros(genome.bases);
			recordBin=63-Long.numberOfLeadingZeros(genome.records);
		}
		final Genome genome;
		final int index, baseBin, recordBin;
		final long rank;
		final byte[] id, fasta;
		Group group;
	}
	private static final class Group{
		Group(final byte[] domain_, final byte[] phylum_){
			domain=domain_;
			phylum=phylum_;
		}
		final byte[] domain, phylum;
		final ArrayList<Member> members=new ArrayList<Member>();
		Member[] sorted;
		int quota;
		long remainder;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Owned state/constants  ----------------*/
	/*--------------------------------------------------------------*/

	private String census, memberships, taxonomy, prefix;
	private int[] sizes;
	private long seed=20261111;
	private Member[] members;
	private final PrintStream outstream;
	private final ArrayList<Group> groups=new ArrayList<Group>();
	private final ObjectMap<Target, Group> byGroup=new ObjectMap<Target, Group>(Target.class, Group.class);
	private final ObjectIntMap<Target> ids=new ObjectIntMap<Target>(Target.class);
	private final Target probe=new Target();
	private final LineParser1 parser=new LineParser1('\t');
	private final ByteBuilder key=new ByteBuilder(), buffer=new ByteBuilder();
	private static final String TAXONOMY_HEADER="input_id\tfasta\tdomain\tphylum\tcorpus_inputs\tcategory\ttaxonomy_status";
	private static final String[] SUFFIXES={".inputs.txt", ".inputs.config", ".flat.config", ".selection.tsv", ".composition.tsv", ".totals.tsv"};
	private static final Comparator<Group> GROUP_ORDER=new Comparator<Group>(){
		public int compare(final Group a, final Group b){
			final int domain=SelectNcrnaBenchmarkInputs.compare(a.domain, b.domain);
			return domain!=0 ? domain : SelectNcrnaBenchmarkInputs.compare(a.phylum, b.phylum);
		}
	};
	private static final Comparator<Group> REMAINDER_ORDER=new Comparator<Group>(){
		public int compare(final Group a, final Group b){
			final int remainder=Long.compare(b.remainder, a.remainder);
			return remainder!=0 ? remainder : GROUP_ORDER.compare(a, b);
		}
	};
	private static final Comparator<Member> MEMBER_ORDER=new Comparator<Member>(){
		public int compare(final Member a, final Member b){
			int value=Integer.compare(a.baseBin, b.baseBin);
			if(value==0){value=Integer.compare(a.recordBin, b.recordBin);}
			if(value==0){value=Long.compare(a.rank, b.rank);}
			return value!=0 ? value : SelectNcrnaBenchmarkInputs.compare(a.id, b.id);
		}
	};
}
