package prok;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.FileFormat;
import map.LongHashSet;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.KillSwitch;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.IntList;
import structures.ListNum;

/** Descriptive context of already graded reference losses. Follows the native
 * A_SampleByteFile parser/I/O lifecycle and the existing ncRNA sequence owner.
 * Reads saved normalized tables and an existing gene-sense reference FASTA;
 * never calls genes, rematches features, or infers a causal failure stage.
 * Classes overlap. Missing context is not a negative observation.
 * @author Diona
 */
public final class ReportNcrnaLossContext{

	/** CLI owns fatal consumer failure so native I/O workers cannot outlive it. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, ReportNcrnaLossContext.class, false);
		try{new ReportNcrnaLossContext(pp.args, pp.outstream).process();}catch(OutOfMemoryError failure){
			KillSwitch.memKill(failure);
			throw failure;
		}catch(RuntimeException | AssertionError failure){
			failure.printStackTrace(pp.outstream);
			KillSwitch.kill("Loss-context report failed; native input workers must not survive a failed consumer");
			throw failure;
		}finally{Shared.closeStream(pp.outstream);}
	}
	/** Finite fixtures deliberately catch invalid-input failures in their own JVM.
	 * This entry does not own a production CLI's fatal worker-lifetime policy. */
	public static void runForTesting(final String[] args){
		final PreParser pp=new PreParser(args, ReportNcrnaLossContext.class, false);
		try{new ReportNcrnaLossContext(pp.args, pp.outstream).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private ReportNcrnaLossContext(final String[] args, final PrintStream stream){
		outstream=stream;
		Shared.setThreads(1);
		final Parser parser=new Parser();
		parser.overwrite=false;
		String orientation=null;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("referencefa")){referenceFa=b;
			}else if(a.equals("contiglengths")){contigLengths=b;
			}else if(a.equals("orientation")){orientation=b;
			}else if(a.equals("expectedrecords")){expectedRecords=Parse.parseKMG(b);
			}else if(a.equals("endpercent")){endPercent=Parse.parseIntKMG(b);
			}else{require(parser.parse(arg, a, b), "Unknown saved loss-context parameter");}
		}
		manifest=parser.in1;
		out=parser.out1;
		require(manifest!=null && out!=null && parser.in2==null && parser.out2==null
			&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Supply the complete explicit case manifest and a fresh output prefix; sampling or overwriting would lose witnesses");
		require(endPercent>0 && endPercent<50, "End-only seed classes need nonoverlapping, declared terminal fractions");
		require(referenceFa==null ? orientation==null && expectedRecords<0 : "sense".equals(orientation) && expectedRecords>0,
			"An existing reference FASTA requires orientation=sense and its accepted expectedrecords; no genome extraction is performed");
		require(Tools.testInputFiles(false, true, manifest, referenceFa, contigLengths), "Declared saved inputs must be readable");
		RrnaResourceIO.fresh(out, ".loci.tsv", ".summary.tsv", ".definitions.tsv");
	}

	private void process(){
		loadCases();
		for(final Case c:cases){
			loadLosses(c);
			loadReferences(c);
			loadQueries(c, 0);
			loadQueries(c, 1);
			classifyNeighbors(c);
		}
		if(contigLengths!=null){loadContigLengths();}
		if(referenceFa!=null && !genes.isEmpty()){loadReferenceSequences();}
		write();
		outstream.println("NCRNA_SAVED_LOSS_CONTEXT_COMPLETE cases="+cases.size()+" unique_reference_sequences="+genes.size()
			+" reference_fasta_records_read="+sequenceRecords+";OVERLAPPING_DESCRIPTIVE_CLASSES;NO_CAUSAL_STAGE_ASSIGNMENT");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Saved tables         ----------------*/
	/*--------------------------------------------------------------*/

	private void loadCases(){
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		scan(manifest, MANIFEST_HEADER, 7, ()->{
			final String family=p.parseString(0), variant=p.parseString(1), seeds=p.parseString(5);
			require(!family.isEmpty() && !variant.isEmpty() && names.add(family+"\t"+variant), "Each family/variant must appear exactly once");
			Family f=families.get(family);
			if(f==null){
				f=new Family(family, seeds);
				for(final String path:seeds.split(",")){
					require(!path.isEmpty(), "Every fixed forward seed path must be explicit");
					final LongHashSet part=ProkObject.loadLongKmers(path, K);
					for(final long word:part.toArray()){f.keys.add(word);}
				}
				require(f.keys.size()>0, "A family needs its unchanged nonempty forward K17 seed set");
				families.put(family, f);
			}else{require(f.seedPaths.equals(seeds), "Variants must name the same fixed family seed resources");}
			final Case c=new Case(f, variant, p.parseString(2), p.parseString(3), p.parseString(4), p.parseInt(6));
			require(c.expectedLosses>=0, "The expected saved loss count must be explicit, including zero");
			cases.add(c);
		});
		require(!cases.isEmpty(), "A context report needs at least one declared family/variant");
		for(final Case c:cases){for(final String prefix:c.prefixes){familyDeclaration(prefix, c.family.name);}}
	}

	private void familyDeclaration(final String prefix, final String family){
		try(Input input=new Input(prefix+".family.tsv")){
			RrnaResourceIO.header(input.reader, "#NCRNA_NORMALIZED_FAMILY_V1", prefix);
			final byte[] row=input.reader.nextLine();
			require(row!=null, "Even an empty query export must declare its biological family");
			p.set(row);
			require(p.terms()==2 && p.termEquals("family", 0) && p.termEquals(family, 1) && input.reader.nextLine()==null,
				"Saved exports must declare exactly the requested biological family");
		}
	}

	private void loadLosses(final Case c){
		scan(c.comparison+".loci.tsv", LOSS_HEADER, 13, ()->{
			if(p.termEquals("GAINED_GRADED_REFERENCE", 6)){return;}
			require(p.termEquals("LOST_GRADED_REFERENCE", 6), "Only the saved comparator's lost/gained classes are understood");
			final Target id=compound(0, 1).copy();
			final int start=p.parseInt(3), stop=p.parseInt(4);
			final char strand=strand(5);
			require(start>0 && stop>=start, "Saved reference coordinates must be one-based inclusive");
			Gene gene=byGene.get(id);
			if(gene==null){
				gene=new Gene(c.family, text(0).copy(), text(1).copy(), text(2).copy(), start, stop, strand);
				byGene.put(id, gene);
				genes.add(gene);
				final Target contig=compound(0, 2);
				GeneGroup group=byContig.get(contig);
				if(group==null){
					group=new GeneGroup();
					byContig.put(contig.copy(), group);
				}
				group.genes.add(gene);
			}else{
				require(gene.family==c.family && gene.seqid.equals(text(2)) && gene.start==start && gene.stop==stop && gene.strand==strand,
					"The same input/reference identity must retain its coordinates, strand and family across variants");
			}
			final Target contig=compound(0, 2);
			Context context=c.contexts.get(contig);
			if(context==null){
				context=new Context();
				c.contexts.put(contig.copy(), context);
			}
			final Loss loss=new Loss(gene, context);
			require(c.byLoss.put(id, loss)==null, "A lost reference must appear exactly once within a comparison");
			c.losses.add(loss);
			context.losses.add(loss);
		});
		require(c.losses.size()==c.expectedLosses, "Saved lost-reference rows must conserve the declared comparator total");
	}

	private void loadReferences(final Case c){
		final IntList lengths=new IntList();
		final ObjectSet<Target> seen=new ObjectSet<Target>(Target.class);
		scan(c.prefixes[0]+".references.tsv", REFERENCE_HEADER, 7, ()->{
			final Target id=compound(0, 1);
			require(seen.add(id.copy()), "Duplicate input-qualified reference IDs would distort lengths and neighborhood counts");
			final Loss loss=c.byLoss.get(id);
			final int start=p.parseInt(3), stop=p.parseInt(4);
			final char strand=strand(5);
			final boolean scored=bool(6);
			require(start>0 && stop>=start, "Reference context intervals must be ordered and positive");
			if(scored){lengths.add(stop-start+1);}
			if(loss!=null){
				final Gene g=loss.gene;
				require(scored && !loss.referenceSeen && g.seqid.equals(text(2)) && start==g.start && stop==g.stop && strand==g.strand,
					"Every lost row must bind its exact original scored reference, not a same-name neighbor");
				loss.referenceSeen=true;
			}
			final Context context=c.contexts.get(compound(0, 2));
			if(context!=null){context.references.add(new Neighbor(text(1).copy(), start, stop, strand, scored));}
		});
		require(lengths.size>0, "Family p1/p99 thresholds require all scored references, not just the lost subset");
		lengths.sort();
		c.scoredReferences=lengths.size;
		c.p1=ReportNcrnaSeedPositions.percentile(lengths.array, lengths.size, 1, 100);
		c.p99=ReportNcrnaSeedPositions.percentile(lengths.array, lengths.size, 99, 100);
		for(final Loss loss:c.losses){require(loss.referenceSeen, "A saved lost reference is missing from the supplied reference view");}
	}

	private void loadQueries(final Case c, final int side){
		assert(side==0 || side==1) : "Only the saved before and after call populations can supply context";
		scan(c.prefixes[side]+".queries.tsv", QUERY_HEADER, 10, ()->{
			require(p.termEquals(c.family.name, 5), "Calls must retain the declared family-union identity");
			final int start=p.parseInt(2), stop=p.parseInt(3);
			final char strand=strand(4);
			require(start>0 && stop>=start, "Saved calls must retain positive inclusive coordinates");
			final boolean tp=p.termEquals("TP", 7);
			require(tp || p.termEquals("MASKED", 7) || p.termEquals("UNMATCHED_DUPLICATE", 7) || p.termEquals("UNMATCHED_OTHER", 7),
				"Call context must use a known normalized outcome");
			final Context context=c.contexts.get(compound(0, 1));
			if(context==null){return;}
			for(final Loss loss:context.losses){
				final Gene g=loss.gene;
				if(tp && g.id.equals(text(9))){
					loss.matches[side]++;
					continue;
				}
				if(start<=g.stop && stop>=g.start){
					loss.overlappingCalls[side]++;
					continue;
				}
				final long gap=gap(g.start, g.stop, start, stop);
				if(gap>g.length()){continue;}
				loss.nearCalls[side]++;
				if(loss.nearestCallGap[side]<0 || gap<loss.nearestCallGap[side]){
					loss.nearestCallGap[side]=gap;
					loss.nearestCallStart[side]=start;
					loss.nearestCallStop[side]=stop;
					loss.nearestCallStrand[side]=strand;
				}
			}
		});
		for(final Loss loss:c.losses){
			require(loss.matches[side]==(side==0 ? 1 : 0), "LOST requires exactly one saved baseline TP and no candidate TP for that reference");
		}
	}

	private void classifyNeighbors(final Case c){
		assert(c.p1>0 && c.p99>=c.p1) : "Length classes use the complete scored family distribution";
		for(final Loss loss:c.losses){
			final Gene g=loss.gene;
			for(final Neighbor other:loss.context.references){
				if(g.id.equals(other.id)){continue;}
				final long gap=gap(g.start, g.stop, other.start, other.stop);
				if(gap>g.length()){continue;}
				loss.nearReferences++;
				if(loss.nearestReference==null || gap<loss.nearestReferenceGap){
					loss.nearestReference=other;
					loss.nearestReferenceGap=gap;
				}
				if(other.scored && other.strand==g.strand && (other.stop<g.start || other.start>g.stop)
					&& gap<=Math.min(g.length(), other.stop-other.start+1)){loss.tandemLike=true;}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Existing evidence       ----------------*/
	/*--------------------------------------------------------------*/

	private void loadContigLengths(){
		scan(contigLengths, LENGTH_HEADER, 3, ()->{
			final int length=p.parseInt(2);
			require(length>0, "A contig length must come from an explicit positive source length, not a maximum feature stop");
			final GeneGroup group=byContig.get(compound(0, 1));
			if(group==null){return;}
			for(final Gene g:group.genes){
				require(g.contigLength<0 && g.stop<=length, "Contig length bindings must be unique and contain the full lost reference");
				g.contigLength=length;
			}
		});
	}

	private void loadReferenceSequences(){
		final LongHashSet distinct=new LongHashSet();
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(
			FileFormat.testInput(referenceFa, FileFormat.FASTA, null, true, false), null, false, -1, false, true, 0))){
			input.stream.start();
			for(ListNum<Read> batch=input.stream.nextList(); batch!=null; batch=input.stream.nextList()){
				for(final Read read:batch){
					sequenceRecords++;
					key.clear();
					decodeAttribute(read.id, "input_id", key);
					key.tab();
					decodeAttribute(read.id, "ID", key);
					probe.set(key.array, 0, key.length());
					final Gene g=byGene.get(probe);
					if(g==null){continue;}
					require(!g.sequenceSeen && read.length()==g.length(), "Saved gene-sense reference sequences must match each target ID exactly once at its full inclusive length");
					g.sequenceSeen=true;
					distinct.clear();
					long word=0;
					int valid=0;
					for(int i=0; i<read.length(); i++){
						final int letter=read.bases[i]&255;
						final int base=letter<AminoAcid.baseToNumber.length ? AminoAcid.baseToNumber[letter] : -1;
						if(base<0){
							word=0;
							valid=0;
							continue;
						}
						word=((word<<2)|base)&KEY_MASK;
						if(++valid<K || !g.family.keys.contains(word)){continue;}
						final int center=i-K/2;
						if(g.minCenter<0){g.minCenter=center;}
						g.maxCenter=center;
						g.keyOccurrences++;
						distinct.add(word);
					}
					g.distinctKeys=distinct.size();
				}
			}
		}
		require(sequenceRecords==expectedRecords, "The saved reference FASTA must conserve its declared complete record population");
		for(final Gene g:genes){require(g.sequenceSeen, "An expected lost reference was not found in the supplied scored-reference FASTA");}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(){
		try(Output loci=new Output(out+".loci.tsv"); Output summary=new Output(out+".summary.tsv"); Output definitions=new Output(out+".definitions.tsv")){
			loci.writer.println(LOCUS_HEADER);
			summary.writer.println("family\tvariant\tclass\tlost_with_class\tassessable_losses\tunique_lost_references\tfraction_all_losses\tqualification");
			definitions.writer.println("key\tvalue");
			definitions.writer.print(row.clear().append("case_manifest\t").append(manifest).nl());
			definitions.writer.print(row.clear().append("reference_fasta\t").append(referenceFa==null ? "NA" : referenceFa).nl());
			definitions.writer.print(row.clear().append("contig_lengths\t").append(contigLengths==null ? "NA" : contigLengths).nl());
			definitions.writer.print(row.clear().append("reference_records_read\t").append(sequenceRecords).nl());
			definitions.writer.println("near_distance\tInter-interval gap at most the lost reference inclusive length; another reference ID or a disjoint call, never the lost locus's own witness");
			definitions.writer.println("tandem_like\tAnother scored, disjoint, same-strand reference within the shorter of the two reference lengths; geometry alone does not prove tandem duplication");
			definitions.writer.println("length_distribution\tNearest-rank p1 and p99 of all scored references in each fixed family view; strict outside comparisons");
			definitions.writer.println("edge\tExact physical contig boundary; a missing right-boundary length is unavailable, not false");
			definitions.writer.print(row.clear().append("seed_end_percent\t").append(endPercent).nl());
			definitions.writer.println("seed_positions\tAll complete fixed forward K17 occurrence centers, zero-based in existing gene-sense reference sequence; repeats retained, ambiguous runs break kmers");
			definitions.writer.println("class_counts\tNonexclusive flags; NONE requires all requested evidence assessed and negative; incomplete negatives are UNASSESSED");
			for(final Case c:cases){
				final long[] counts=new long[CLASSES.length], assessed=new long[CLASSES.length];
				for(final Loss loss:c.losses){
					final Gene g=loss.gene;
					final int[] flags=flags(loss, c);
					for(int i=0; i<flags.length; i++){
						if(flags[i]>=0){assessed[i]++;}
						if(flags[i]>0){counts[i]++;}
					}
					row.clear().append(c.family.name).tab().append(c.variant).tab();
					append(row, g.input).tab();
					append(row, g.id).tab();
					append(row, g.seqid).tab();
					row.append(g.start).tab().append(g.stop).tab().append(g.strand).tab().append(g.length()).tab().append(c.scoredReferences).tab().append(c.p1).tab().append(c.p99).tab();
					optional(row, g.contigLength).tab().append(g.start-1).tab();
					optional(row, g.contigLength<0 ? -1 : g.contigLength-g.stop).tab();
					row.append(loss.nearReferences).tab();
					optional(row, loss.nearestReferenceGap).tab();
					if(loss.nearestReference==null){row.append("NA");}else{append(row, loss.nearestReference.id);}
					for(int side=0; side<2; side++){
						row.tab().append(loss.nearCalls[side]).tab();
						optional(row, loss.nearestCallGap[side]).tab();
						optional(row, loss.nearestCallStart[side]).tab();
						optional(row, loss.nearestCallStop[side]).tab().append(loss.overlappingCalls[side]);
					}
					row.tab();
					optional(row, g.sequenceSeen ? g.keyOccurrences : -1).tab();
					optional(row, g.sequenceSeen ? g.distinctKeys : -1).tab();
					optional(row, g.minCenter).tab();
					optional(row, g.maxCenter);
					for(final int flag:flags){row.tab().append(flag<0 ? "NA" : flag>0 ? "true" : "false");}
					loci.writer.print(row.append("\t").append(QUALIFICATION).nl());
				}
				for(int i=0; i<CLASSES.length; i++){
					row.clear().append(c.family.name).tab().append(c.variant).tab().append(CLASSES[i]).tab().append(counts[i]).tab().append(assessed[i]).tab().append(c.losses.size()).tab();
					if(c.losses.isEmpty()){row.append("NA");}else{row.append(counts[i]/(double)c.losses.size(), 8);}
					summary.writer.print(row.tab().append(QUALIFICATION).nl());
				}
			}
		}
	}

	private int[] flags(final Loss loss, final Case c){
		final Gene g=loss.gene;
		final int edge=g.start==1 || g.contigLength==g.stop ? 1 : g.contigLength<0 ? -1 : 0;
		final boolean hasKeys=g.keyOccurrences>0;
		final int[] result={loss.tandemLike ? 1 : 0, loss.nearReferences>0 ? 1 : 0,
			loss.nearCalls[0]>0 ? 1 : 0, loss.nearCalls[1]>0 ? 1 : 0, edge,
			g.length()<c.p1 ? 1 : 0, g.length()>c.p99 ? 1 : 0,
			!g.sequenceSeen ? -1 : hasKeys && (long)g.maxCenter*100<=(long)endPercent*(g.length()-1) ? 1 : 0,
			!g.sequenceSeen ? -1 : hasKeys && (long)g.minCenter*100>=(100L-endPercent)*(g.length()-1) ? 1 : 0,
			!g.sequenceSeen ? -1 : hasKeys ? 0 : 1, 0, 0};
		boolean any=false, missing=false;
		for(int i=0; i<result.length-2; i++){
			any|=result[i]>0;
			missing|=result[i]<0;
		}
		result[10]=!any && !missing ? 1 : 0;
		result[11]=!any && missing ? 1 : 0;
		assert(result.length==CLASSES.length && result[10]+result[11]<=1) : "NONE and unassessed cannot double-count a negative lost-reference context";
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Shared parsing        ----------------*/
	/*--------------------------------------------------------------*/

	private void scan(final String path, final String header, final int fields, final Row action){
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, header, path);
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				p.set(line);
				if(p.terms()!=fields){throw new IllegalArgumentException("Saved table schema changed: "+path);}
				action.read();
			}
		}
	}
	private Target text(final int field){
		probe.set(p, field);
		return probe;
	}
	private Target compound(final int first, final int second){
		key.clear();
		p.setBounds(first);
		require(p.b()>p.a(), "Input identities cannot be empty");
		key.append(p.line(), p.a(), p.b()-p.a()).tab();
		p.setBounds(second);
		require(p.b()>p.a(), "Reference/contig identities cannot be empty");
		key.append(p.line(), p.a(), p.b()-p.a());
		probe.set(key.array, 0, key.length());
		return probe;
	}
	private char strand(final int field){
		require(p.termEquals("+", field) || p.termEquals("-", field), "Reference and call strand must be explicit");
		return p.termEquals("+", field) ? '+' : '-';
	}
	private boolean bool(final int field){
		require(p.termEquals("true", field) || p.termEquals("false", field), "Scored disposition must be explicit true/false");
		return p.termEquals("true", field);
	}
	/** Inclusive intervals separated by no bases have zero gap; overlap is also zero. */
	public static long gap(final int a, final int z, final int b, final int y){
		assert(a>0 && z>=a && b>0 && y>=b) : "Context gaps operate on validated one-based inclusive intervals";
		return b>z ? b-(long)z-1 : a>y ? a-(long)y-1 : 0;
	}
	/** Decode each attribute exactly once. Literal '+' stays '+', unlike form URLs. */
	public static void decodeAttribute(final String title, final String name, final ByteBuilder output){
		assert(title!=null && name!=null && output!=null) : "Saved gene headers and the reusable destination must exist before ID binding";
		final String marker="ID".equals(name) ? ID_MARKER : "input_id".equals(name) ? INPUT_MARKER : null;
		require(marker!=null, "Only fixed reference ID and input_id attributes are used for sequence binding");
		int start=title.indexOf(marker);
		while(start>0 && title.charAt(start-1)!=';'){start=title.indexOf(marker, start+1);}
		if(start<0){throw new IllegalArgumentException("Saved reference FASTA lacks required attribute: "+name);}
		start+=marker.length();
		int end=start;
		while(end<title.length() && title.charAt(end)!=';' && title.charAt(end)>32){end++;}
		require(end>start, "Saved reference FASTA identity attributes cannot be empty");
		for(int i=start; i<end; i++){
			final char c=title.charAt(i);
			if(c=='%'){
				require(i+2<end, "Truncated percent escape in reference identity");
				final int a=Character.digit(title.charAt(++i), 16), b=Character.digit(title.charAt(++i), 16);
				require(a>=0 && b>=0, "Invalid percent escape in reference identity");
				output.append((byte)((a<<4)|b));
			}else{
				require(c<128, "Non-ASCII GFF identity bytes must be percent-encoded");
				output.append((byte)c);
			}
		}
	}
	private static ByteBuilder append(final ByteBuilder output, final Target value){return output.append(value.bytes, value.from, value.to-value.from);}
	private static ByteBuilder optional(final ByteBuilder output, final long value){return value<0 ? output.append("NA") : output.append(value);}
	private static void require(final boolean valid, final String message){if(!valid){throw new IllegalArgumentException(message);}}
	private interface Row{void read();}

	/*--------------------------------------------------------------*/
	/*----------------         Retained state       ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Family{
		Family(final String name_, final String paths){
			name=name_;
			seedPaths=paths;
		}
		final String name, seedPaths;
		final LongHashSet keys=new LongHashSet();
	}
	private static final class Case{
		Case(final Family family_, final String variant_, final String before, final String after, final String comparison_, final int expected){
			family=family_;
			variant=variant_;
			prefixes=new String[]{before, after};
			comparison=comparison_;
			expectedLosses=expected;
		}
		final Family family;
		final String variant, comparison;
		final String[] prefixes;
		final int expectedLosses;
		final ArrayList<Loss> losses=new ArrayList<Loss>();
		final ObjectMap<Target,Loss> byLoss=new ObjectMap<Target,Loss>(Target.class, Loss.class);
		final ObjectMap<Target,Context> contexts=new ObjectMap<Target,Context>(Target.class, Context.class);
		int scoredReferences, p1, p99;
	}
	private static final class Gene{
		Gene(final Family family_, final Target input_, final Target id_, final Target seqid_, final int start_, final int stop_, final char strand_){
			family=family_;
			input=input_;
			id=id_;
			seqid=seqid_;
			start=start_;
			stop=stop_;
			strand=strand_;
		}
		int length(){return stop-start+1;}
		final Family family;
		final Target input, id, seqid;
		final int start, stop;
		final char strand;
		int contigLength=-1, minCenter=-1, maxCenter=-1, distinctKeys;
		long keyOccurrences;
		boolean sequenceSeen;
	}
	private static final class GeneGroup{final ArrayList<Gene> genes=new ArrayList<Gene>();}
	private static final class Context{
		final ArrayList<Loss> losses=new ArrayList<Loss>();
		final ArrayList<Neighbor> references=new ArrayList<Neighbor>();
	}
	private static final class Neighbor{
		Neighbor(final Target id_, final int start_, final int stop_, final char strand_, final boolean scored_){
			id=id_;
			start=start_;
			stop=stop_;
			strand=strand_;
			scored=scored_;
		}
		final Target id;
		final int start, stop;
		final char strand;
		final boolean scored;
	}
	private static final class Loss{
		Loss(final Gene gene_, final Context context_){
			gene=gene_;
			context=context_;
		}
		final Gene gene;
		final Context context;
		final int[] matches=new int[2], nearCalls=new int[2], overlappingCalls=new int[2];
		final int[] nearestCallStart={-1,-1}, nearestCallStop={-1,-1};
		final char[] nearestCallStrand=new char[2];
		final long[] nearestCallGap={-1,-1};
		Neighbor nearestReference;
		long nearestReferenceGap=-1;
		int nearReferences;
		boolean referenceSeen, tandemLike;
	}

	private String manifest, referenceFa, contigLengths, out;
	private long expectedRecords=-1, sequenceRecords;
	private int endPercent=25;
	private final PrintStream outstream;
	private final ArrayList<Case> cases=new ArrayList<Case>();
	private final ArrayList<Gene> genes=new ArrayList<Gene>();
	private final ObjectMap<String,Family> families=new ObjectMap<String,Family>(String.class, Family.class);
	private final ObjectMap<Target,Gene> byGene=new ObjectMap<Target,Gene>(Target.class, Gene.class);
	private final ObjectMap<Target,GeneGroup> byContig=new ObjectMap<Target,GeneGroup>(Target.class, GeneGroup.class);
	private final LineParser1 p=new LineParser1('\t');
	private final Target probe=new Target();
	private final ByteBuilder key=new ByteBuilder(), row=new ByteBuilder();
	private static final int K=17;
	private static final long KEY_MASK=(1L<<(2*K))-1;
	private static final String ID_MARKER="ID=", INPUT_MARKER="input_id=";
	public static final String MANIFEST_HEADER="family\tvariant\tbefore\tafter\tcomparison\tseeds\texpected_losses";
	public static final String REFERENCE_HEADER="input_id\treference_id\tseqid\tstart\tstop\tstrand\tscored";
	public static final String QUERY_HEADER="input_id\tseqid\tstart\tstop\tstrand\tfamily_union\tmodel\toutcome\tscored_eligible\tmatched_reference_id";
	public static final String LOSS_HEADER="input_id\treference_id\tseqid\tstart\tstop\tstrand\tchange\tpresence_class\tbaseline_family_union\tbaseline_model\tbaseline_model_in_subset\tcandidate_family_union\tcandidate_model";
	public static final String LENGTH_HEADER="input_id\tseqid\tlength";
	private static final String QUALIFICATION="SAVED_REFERENCE_CONTEXT;OVERLAPPING_FLAGS;GEOMETRY_NOT_CAUSAL_FAILURE_STAGE";
	private static final String[] CLASSES={"TANDEM_LIKE_SCORED_REFERENCE", "NEAR_OTHER_REFERENCE", "NEAR_DISJOINT_CALL_BEFORE", "NEAR_DISJOINT_CALL_AFTER",
		"AT_CONTIG_EDGE", "SHORT_BELOW_P1", "LONG_ABOVE_P99", "KEYS_5PRIME_ONLY", "KEYS_3PRIME_ONLY", "NO_FORWARD_SEED", "NONE_OF_THESE", "UNASSESSED_WITHOUT_POSITIVE_CLASS"};
	private static final String LOCUS_HEADER="family\tvariant\tinput_id\treference_id\tseqid\tstart\tstop\tstrand\treference_bp\tscored_length_population\tlength_p1\tlength_p99\tcontig_length\tleft_edge_distance\tright_edge_distance\tnear_reference_count\tnearest_reference_gap\tnearest_reference_id\tnear_call_before_count\tnearest_call_before_gap\tnearest_call_before_start\tnearest_call_before_stop\toverlapping_other_calls_before\tnear_call_after_count\tnearest_call_after_gap\tnearest_call_after_start\tnearest_call_after_stop\toverlapping_other_calls_after\tseed_occurrences\tdistinct_forward_seeds\tseed_min_center0_sense\tseed_max_center0_sense\ttandem_like_scored_reference\tnear_other_reference\tnear_disjoint_call_before\tnear_disjoint_call_after\tat_contig_edge\tshort_below_p1\tlong_above_p99\tkeys_5prime_only\tkeys_3prime_only\tno_forward_seed\tnone_of_these\tunassessed_without_positive_class\tqualification";
}
