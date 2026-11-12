package prok;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaJointPolicy.Member;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaSweep.Table;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Reports recall/work frontiers from completed native search metadata only.
 * No caller, alignment, genome, GFF or taxonomy-tree access occurs here. Groups
 * must share the exact population, reference, runtime declarations and caller
 * settings. Work counters are logical proxies, not whole-caller time. Retained
 * per-input vectors prove only equality of the reported counters, not loci.
 * One instance owns all readers, writers and retained state on one thread.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSweepFrontier{

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final NcrnaSweepFrontier report=new NcrnaSweepFrontier(args);
		try{
			report.process();
			timer.stop();
			report.outstream.println("Reported saved policies; no caller work; time: "+timer);
		}finally{Shared.closeStream(report.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Adapts template.A_SampleByteFile's parser and native table I/O ownership. */
	public NcrnaSweepFrontier(final String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("taxonomy")){taxonomy=b;}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		manifest=parser.in1;
		prefix=parser.out1;
		require(manifest!=null && taxonomy!=null && prefix!=null && !parser.append && parser.in2==null && parser.out2==null,
			"Supply in=search_manifest, taxonomy=full_selection_ledger, out=prefix; no append/paired input");
		require(Tools.testInputFiles(false, true, manifest, taxonomy), "Cannot read saved metadata declarations");
		final String[] paths=new String[SUFFIXES.length+2];
		paths[0]=manifest;
		paths[1]=taxonomy;
		for(int i=0; i<SUFFIXES.length; i++){paths[i+2]=prefix+SUFFIXES[i];}
		require(Tools.testForDuplicateFiles(true, paths), "Report inputs and outputs must be distinct");
		outputs=new FileFormat[SUFFIXES.length];
		for(int i=0; i<outputs.length; i++){
			require(Tools.testOutputFiles(parser.overwrite, false, false, paths[i+2]), "Cannot write report output: "+paths[i+2]);
			outputs[i]=FileFormat.testOutput(paths[i+2], FileFormat.TXT, null, true, parser.overwrite, false, false);
		}
	}

	private void process(){
		readTaxonomy();
		final ArrayList<Group> groups=readManifest();
		try(Output all=new Output(outputs[0]); Output front=new Output(outputs[1]); Output strata=new Output(outputs[2]);
			Output taxa=new Output(outputs[3]); Output groupOut=new Output(outputs[4]); Output gaps=new Output(outputs[5])){
			all.writer.println(BREADTH_NOTE);
			all.writer.println(HEADER);
			front.writer.println(BREADTH_NOTE);
			front.writer.println(HEADER);
			strata.writer.println("group\tfamily\tsearch\tround\tpolicy\tcategory\tfirst_gain_inputs\tlast_loss_inputs");
			taxa.writer.println("input_id\tfasta\tdomain\tphylum\tcorpus_inputs\tcategory\ttaxonomy_status");
			groupOut.writer.println("group\tfamily\tsearch\troot\tidentity\tcaller\treference\tfamilymap\tsettings");
			gaps.writer.println("group\tsearch\tround\tstatus\tpath");
			writeTaxonomy(taxa);
			for(final Group group:groups){
				for(final Search search:group.searches){
					group.bind(search);
					writeGroup(groupOut, group, search);
					readSearch(group, search, gaps);
				}
				require(!group.points.isEmpty(), "No complete saved policy rows in group "+group.name);
				for(final Point point:group.points){
					if(point.duplicate==null){
						point.tpFront=point.values[BASES]>=0;
						point.breadthFront=point.values[BASES]>=0;
						for(final Point other:group.points){
							if(other.duplicate!=null){continue;}
							point.tpFront&=!dominates(other.values[TP], other.values[BASES], other.first(), other.last(),
								point.values[TP], point.values[BASES], point.first(), point.last(), false);
							point.breadthFront&=!dominates(other.values[TP], other.values[BASES], other.first(), other.last(),
								point.values[TP], point.values[BASES], point.first(), point.last(), true);
						}
					}
				}
				for(final Point point:group.points){
					writePoint(all, group, point);
					if(point.duplicate==null){writeStrata(strata, group, point);}
					if(point.duplicate==null && (point.tpFront || point.breadthFront || point.baseline)){
						writePoint(front, group, point);
					}
				}
				// Group-local per-input vectors are released before the next population.
				group.points.clear();
				group.byTuple.clear();
				group.baseline=null;
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Corpus Taxonomy       ----------------*/
	/*--------------------------------------------------------------*/

	/** Original input IDs are unique; membership rows must not inflate frequency. */
	private void readTaxonomy(){
		try(Table table=new Table(taxonomy)){
			final int id=table.column("input_id"), fasta=table.column("fasta");
			final int domain=table.column("domain"), phylum=table.column("phylum"), status=table.column("taxonomy_status");
			while(table.next()){
				final String name=table.text(id), path=table.text(fasta), d=table.literal(domain), p=table.literal(phylum);
				require(!name.isEmpty() && !byInput.contains(name), "Full-corpus taxonomy requires distinct original input IDs: "+name);
				final String key=d+"/"+p;
				Lineage lineage=byLineage.get(key);
				if(lineage==null){
					final int category=unknown(d) || unknown(p) ? 2
						: p.equalsIgnoreCase("UNCLASSIFIED") || d.equalsIgnoreCase("UNCLASSIFIED") ? 3 : 1;
					lineage=new Lineage(d, p, category);
					byLineage.put(key, lineage);
					lineages.add(lineage);
				}
				lineage.frequency++;
				final Taxon taxon=new Taxon(name, path, table.literal(status), lineage);
				byInput.put(name, taxon);
				taxa.add(taxon);
			}
		}
		require(!taxa.isEmpty(), "Corpus taxonomy cannot be empty");
		final long[] frequencies=new long[lineages.size()];
		int informative=0;
		for(final Lineage lineage:lineages){if(lineage.category==1){frequencies[informative++]=lineage.frequency;}}
		if(informative>0){
			Arrays.sort(frequencies, 0, informative);
			// Nearest-rank lower quartile over lineage labels, including all ties.
			final long cutoff=frequencies[(informative-1)/4];
			for(final Lineage lineage:lineages){
				if(lineage.category==1 && lineage.frequency<=cutoff){lineage.category=0;}
			}
		}
	}

	private void writeTaxonomy(final Output out){
		for(final Taxon taxon:taxa){
			field(buffer.clear(), taxon.id).tab();
			field(buffer, taxon.fasta).tab();
			field(buffer, taxon.lineage.domain).tab();
			field(buffer, taxon.lineage.phylum).tab().append(taxon.lineage.frequency).tab()
				.append(CATEGORIES[taxon.lineage.category]).tab();
			out.writer.print(field(buffer, taxon.status).nl());
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Saved Searches       ----------------*/
	/*--------------------------------------------------------------*/

	private ArrayList<Group> readManifest(){
		final ArrayList<Group> groups=new ArrayList<Group>();
		final ObjectMap<String, Group> byGroup=new ObjectMap<String, Group>(String.class, Group.class);
		final ObjectIntMap<String> roots=new ObjectIntMap<String>(String.class);
		try(Table table=new Table(manifest)){
			final int g=table.column("group"), f=table.column("family"), s=table.column("search"), r=table.column("root");
			while(table.next()){
				final String name=table.literal(g), family=table.literal(f), label=table.literal(s), root=table.text(r);
				require(!name.isEmpty() && !family.isEmpty() && !label.isEmpty() && !roots.contains(root), "Search roots must be unique and labeled");
				roots.put(root, roots.size());
				Group group=byGroup.get(name);
				if(group==null){
					group=new Group(name, family);
					groups.add(group);
					byGroup.put(name, group);
				}
				require(group.family.equals(family), "A comparable group cannot mix biological families");
				group.searches.add(new Search(label, new File(root)));
			}
		}
		require(!groups.isEmpty(), "Supply at least one saved search");
		return groups;
	}

	/** Round directories are numeric and contiguous; no recursive scan occurs. */
	private void readSearch(final Group group, final Search search, final Output gaps){
		final File[] children=search.root.listFiles();
		require(children!=null, "Cannot list declared search root: "+search.root);
		int last=-1;
		for(final File child:children){
			if(child.isDirectory() && child.getName().matches("round[0-9]+")){
				last=Math.max(last, Integer.parseInt(child.getName().substring(5)));
			}
		}
		require(last>=0, "No saved round directories at "+search.root);
		for(int round=0; round<=last; round++){
			final File directory=new File(search.root, "round"+round);
			if(!new File(directory, "complete.txt").isFile()){
				field(buffer.clear(), group.name).tab();
				field(buffer, search.name).tab().append(round).append("\tINCOMPLETE_ROUND_NOT_REPLAYED\t");
				gaps.writer.print(field(buffer, directory.getPath()).nl());
				continue;
			}
			final String marker=text(new File(directory, "complete.txt"));
			require(marker.equals("native\tEVALUATION_RETURNED_AND_SUMMARY_PARSED\n")
				|| marker.equals("native\tEVALUATION_RETURNED_AND_COVERAGE_VALIDATED\n"), "Only completed native rounds belong in the biological report");
			final NcrnaJointPolicy[] policies=NcrnaJointPolicy.load(new File(directory, "policies.tsv").getPath());
			require(policies[0].name.equals("shipped"), "Every saved round must begin with shipped baseline");
			final Point[] points=new Point[policies.length];
			for(int p=0; p<points.length; p++){points[p]=new Point(search.name, round, policies[p], group.bindings.length);}
			readInputs(group, directory, points);
			readSummary(directory, points);
			for(final Point point:points){
				if(group.baseline==null){group.baseline=point;}
				if(point.baseline){
					require(point.tuple.equals(group.baseline.tuple) && point.same(group.baseline),
						"Shipped counters changed in comparable group; split or investigate without rerunning: "+group.name);
				}
				final PointList previous=group.byTuple.get(point.tuple);
				if(previous==null){
					final PointList list=new PointList();
					list.add(point);
					group.byTuple.put(point.tuple, list);
				}else{
					for(final Point candidate:previous){
						if(point.same(candidate)){
							point.duplicate=candidate;
							break;
						}
					}
					if(point.duplicate==null){
						point.conflict=true;
						for(final Point candidate:previous){candidate.conflict=true;}
						previous.add(point);
					}else{point.inputValues=null;}
				}
				group.points.add(point);
			}
		}
	}

	/** Native input order and all reported vector fields are checked once per round. */
	private void readInputs(final Group group, final File directory, final Point[] points){
		try(Table table=new Table(new File(directory, "sweep.inputs.tsv").getPath())){
			final int id=table.column("input_id"), policy=table.column("policy");
			final int[] columns=columns(table);
			for(int i=0; i<group.bindings.length; i++){
				final Binding binding=group.bindings[i];
				final Taxon taxon=byInput.get(binding.id);
				require(taxon!=null && taxon.fasta.equals(binding.fasta), "Reference input must join exact original taxonomy ID and FASTA: "+binding.id);
				for(int p=0; p<points.length; p++){
					final Point point=points[p];
					require(table.next() && table.matches(group.tokens[i], id) && table.matches(point.token, policy),
						"Missing, repeated or unexpected native input/policy row in "+directory);
					final int offset=i*METRICS.length;
					for(int c=0; c<columns.length; c++){
						final long value=number(table, columns[c]);
						require(c>GAINED || value>=0, "Detection counts cannot be unmeasured: "+METRICS[c]);
						point.inputValues[offset+c]=value;
						point.values[c]=sum(point.values[c], value);
					}
					final long tp=point.inputValues[offset+TP], scored=point.inputValues[offset+SCORED];
					final long baseline=points[0].inputValues[offset+TP];
					require(Math.addExact(tp, point.inputValues[offset+FN])==scored
						&& scored==points[0].inputValues[offset+SCORED], "Native input TP+FN and fixed scored denominator must conserve");
					final long lost=point.inputValues[offset+LOST], gained=point.inputValues[offset+GAINED];
					require(lost<=baseline && gained<=scored-baseline && Math.addExact(baseline-lost, gained)==tp,
						"Saved per-input paired changes must explain policy TP");
					if(p==0){require(lost==0 && gained==0, "Shipped input must be its own change baseline");}
					if(baseline==0 && tp>0){point.first[taxon.lineage.category]++;}
					if(baseline>0 && tp==0){point.last[taxon.lineage.category]++;}
				}
			}
			require(!table.next(), "Extra input/policy rows cannot be ignored");
		}
	}

	private static void readSummary(final File directory, final Point[] points){
		try(Table table=new Table(new File(directory, "sweep.summary.tsv").getPath())){
			final int policy=table.column("policy");
			final int[] columns=columns(table);
			for(final Point point:points){
				require(table.next() && table.matches(point.token, policy), "Every saved policy needs exactly one native summary row");
				for(int c=0; c<columns.length; c++){
					require(number(table, columns[c])==point.values[c], "Native per-input totals disagree with summary: "+METRICS[c]);
				}
			}
			require(!table.next(), "Unexpected extra summary policy");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Frontier Output       ----------------*/
	/*--------------------------------------------------------------*/

	/** Strict Pareto dominance; FP/CDS are reports and never selection axes. */
	public static boolean dominates(final long tpA, final long basesA, final long firstA, final long lastA,
			final long tpB, final long basesB, final long firstB, final long lastB, final boolean breadth){
		require(tpA>=0 && tpB>=0 && firstA>=0 && firstB>=0 && lastA>=0 && lastB>=0,
			"Frontier detection and input counts must be measured nonnegative values");
		if(basesA<0 || basesB<0){return false;}
		return tpA>=tpB && basesA<=basesB && (!breadth || firstA>=firstB && lastA<=lastB)
			&& (tpA>tpB || basesA<basesB || breadth && (firstA>firstB || lastA<lastB));
	}

	private void writePoint(final Output out, final Group group, final Point point){
		field(buffer.clear(), group.name).tab();
		field(buffer, group.family).tab();
		field(buffer, point.search).tab().append(point.round).tab();
		field(buffer, point.policy.name).tab().append(point.baseline).tab().append(group.bindings.length);
		for(final long value:point.values){
			buffer.tab();
			count(buffer, value);
		}
		buffer.tab();
		percent(buffer, point.values[TP], point.values[SCORED], false);
		buffer.tab();
		percent(buffer, point.values[TP], Math.addExact(point.values[TP], point.values[FP]), false);
		buffer.tab().append(point.first()).tab().append(point.last()).tab();
		percent(buffer, point.values[ALIGN], group.baseline.values[ALIGN], true);
		buffer.tab();
		percent(buffer, point.values[BASES], group.baseline.values[BASES], true);
		buffer.tab().append(point.values[BASES]<0 ? "NA" : Boolean.toString(point.duplicate==null ? point.tpFront : point.duplicate.tpFront));
		buffer.tab().append(point.values[BASES]<0 ? "NA" : Boolean.toString(point.duplicate==null ? point.breadthFront : point.duplicate.breadthFront));
		buffer.tab().append(point.conflict || point.duplicate!=null && point.duplicate.conflict ? "CONFLICTING_TUPLE"
			: point.duplicate!=null ? "CONSISTENT_DUPLICATE" : "UNIQUE_MEASUREMENT").tab();
		field(buffer, point.duplicate==null ? "." : point.duplicate.location()).tab();
		field(buffer, point.tuple).append("\tUNMEASURED\n");
		out.writer.print(buffer);
	}

	private void writeStrata(final Output out, final Group group, final Point point){
		for(int c=0; c<CATEGORIES.length; c++){
			field(buffer.clear(), group.name).tab();
			field(buffer, group.family).tab();
			field(buffer, point.search).tab().append(point.round).tab();
			field(buffer, point.policy.name).tab().append(CATEGORIES[c]).tab().append(point.first[c]).tab().append(point.last[c]);
			out.writer.print(buffer.nl());
		}
	}

	private void writeGroup(final Output out, final Group group, final Search search){
		field(buffer.clear(), group.name).tab();
		field(buffer, group.family).tab();
		field(buffer, search.name).tab();
		field(buffer, search.root.getPath());
		for(final String name:new String[]{"identity.tsv", "caller.config", "reference.tsv", "familymap.tsv"}){
			buffer.tab();
			field(buffer, new File(search.root, name).getPath());
		}
		out.writer.print(field(buffer.tab(), group.settings).nl());
	}

	/** Percentages use measured longs; unavailable/zero denominators remain NA. */
	private static void percent(final ByteBuilder out, final long value, final long denominator, final boolean change){
		if(value<0 || denominator<=0){out.append("NA");}
		else{out.append(100*(value/(double)denominator-(change ? 1 : 0)), 6);}
	}
	private static void count(final ByteBuilder out, final long value){
		if(value<0){out.append("NA");}
		else{out.append(value);}
	}
	private static long sum(final long a, final long b){return a<0 || b<0 ? -1 : Math.addExact(a, b);}
	private static boolean unknown(final String name){
		return name.isEmpty() || name.equals(".") || name.equalsIgnoreCase("UNKNOWN") || name.equalsIgnoreCase("NA");
	}
	private static long number(final Table table, final int column){
		return column<0 || table.parser.termEquals("NA", column) ? -1 : NcrnaSearchCoverage.count(table, column);
	}
	private static int[] columns(final Table table){
		assert(table!=null) : "Saved frontier columns require an owned header-indexed table";
		final int[] columns=new int[METRICS.length];
		for(int c=0; c<columns.length; c++){
			columns[c]=table.columns.contains(METRICS[c]) ? table.column(METRICS[c]) : -1;
			require(c>GAINED || columns[c]>=0, "Detection fields must exist in saved native metadata: "+METRICS[c]);
		}
		return columns;
	}
	private static String tuple(final NcrnaJointPolicy policy){
		assert(policy!=null && policy.size()>0) : "A complete joint policy supplies every member tuple used for deduplication";
		final StringBuilder key=new StringBuilder();
		for(int m=0; m<policy.size(); m++){
			final Member member=policy.member(m);
			key.append(member.name).append(':').append(member.adaptive).append(':').append(member.topN).append(':')
				.append(member.floor).append(':').append(member.minHits).append(':').append(member.topFrac).append(':')
				.append(member.qFrac).append(':').append(member.patience).append(';');
		}
		return key.toString();
	}
	private static String text(final File file){
		final ByteBuilder text=new ByteBuilder();
		try(Input in=new Input(file.getPath())){
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){text.append(line).nl();}
		}
		return new String(text.toBytes(), StandardCharsets.UTF_8);
	}
	private static ByteBuilder field(final ByteBuilder out, final String text){return BuildNcrnaReference.field(out, text);}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Row Owners          ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Group{
		Group(final String name_, final String family_){
			name=name_;
			family=family_;
		}
		void bind(final Search search){
			final String[] names={"identity.tsv", "caller.config", "reference.tsv", "familymap.tsv"};
			final String[] found=new String[names.length];
			for(int i=0; i<names.length; i++){found[i]=text(new File(search.root, names[i]));}
			final StringBuilder controls=new StringBuilder();
			final ObjectMap<String, String> contract=new ObjectMap<String, String>(String.class, String.class);
			for(final String line:text(new File(search.root, "contract.txt")).split("\n")){
				final int equal=line.indexOf('=');
				require(equal>0, "Search contract requires key=value declarations");
				final String key=line.substring(0, equal);
				require(!contract.contains(key), "Repeated search contract key: "+key);
				contract.put(key, line.substring(equal+1));
			}
			for(final String key:COMPARABLE_SETTINGS){
				require(contract.contains(key), "Missing caller-comparability setting: "+key);
				controls.append(key).append('=').append(contract.get(key)).append(';');
			}
			require("native".equals(contract.get("backend")), "Synthetic evaluators cannot become biological measurements");
			if(signature==null){
				signature=found;
				settings=controls.toString();
				final String[] caller=Parser.parseConfig(new String[]{"config="
					+NcrnaCallerProcess.configurationPath(new File(search.root, "caller.config"))});
				bindings=NcrnaSearchCoverage.bindInputs(caller, new File(search.root, "reference.tsv").getPath(), true);
				tokens=new byte[bindings.length][];
				for(int i=0; i<tokens.length; i++){tokens[i]=NcrnaSearchCoverage.token(bindings[i].id);}
			}else{
				require(Arrays.equals(signature, found) && settings.equals(controls.toString()),
					"Do not compare different panels/references/runtimes/caller settings inside group "+name);
			}
		}
		final String name, family;
		String[] signature;
		String settings;
		Binding[] bindings;
		byte[][] tokens;
		Point baseline;
		final ArrayList<Search> searches=new ArrayList<Search>();
		final ArrayList<Point> points=new ArrayList<Point>();
		final ObjectMap<String, PointList> byTuple=new ObjectMap<String, PointList>(String.class, PointList.class);
	}
	/** One retained list per distinct policy tuple, including conflicting outcomes. */
	private static final class PointList extends ArrayList<Point>{
		private static final long serialVersionUID=1L;
	}
	private static final class Search{
		Search(final String name_, final File root_){
			name=name_;
			root=root_;
		}
		final String name;
		final File root;
	}
	private static final class Point{
		Point(final String search_, final int round_, final NcrnaJointPolicy policy_, final int inputs){
			search=search_;
			round=round_;
			policy=policy_;
			token=NcrnaSearchCoverage.token(policy.name);
			tuple=tuple(policy);
			baseline=policy.name.equals("shipped");
			inputValues=new long[Math.multiplyExact(inputs, METRICS.length)];
		}
		boolean same(final Point other){
			return Arrays.equals(values, other.values) && Arrays.equals(inputValues, other.inputValues)
				&& Arrays.equals(first, other.first) && Arrays.equals(last, other.last);
		}
		long first(){return first[0]+first[1]+first[2]+first[3];}
		long last(){return last[0]+last[1]+last[2]+last[3];}
		String location(){return search+"/round"+round+"/"+policy.name;}
		final String search, tuple;
		final int round;
		final NcrnaJointPolicy policy;
		final byte[] token;
		final boolean baseline;
		final long[] values=new long[METRICS.length], first=new long[4], last=new long[4];
		long[] inputValues;
		Point duplicate;
		boolean conflict, tpFront, breadthFront;
	}
	private static final class Lineage{
		Lineage(final String domain_, final String phylum_, final int category_){
			domain=domain_;
			phylum=phylum_;
			category=category_;
		}
		final String domain, phylum;
		long frequency;
		int category;
	}
	private static final class Taxon{
		Taxon(final String id_, final String fasta_, final String status_, final Lineage lineage_){
			id=id_;
			fasta=fasta_;
			status=status_;
			lineage=lineage_;
		}
		final String id, fasta, status;
		final Lineage lineage;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String manifest, taxonomy, prefix;
	private final PrintStream outstream;
	private final FileFormat[] outputs;
	private final ByteBuilder buffer=new ByteBuilder();
	private final ObjectMap<String, Taxon> byInput=new ObjectMap<String, Taxon>(String.class, Taxon.class);
	private final ObjectMap<String, Lineage> byLineage=new ObjectMap<String, Lineage>(String.class, Lineage.class);
	private final ArrayList<Taxon> taxa=new ArrayList<Taxon>();
	private final ArrayList<Lineage> lineages=new ArrayList<Lineage>();
	private static final int TP=0, FP=1, FN=2, SCORED=3, LOST=4, GAINED=5, ALIGN=6, BASES=7;
	private static final String[] METRICS={"TP", "FP", "FN", "scored_references", "lost_baseline_TP", "gained_TP",
		"kernel_logical_alignments", "aligned_model_bases", "CDS"};
	private static final String[] COMPARABLE_SETTINGS={"schema", "patience", "backend", "target", "grading", "cachewindows",
		"cachelength", "reusecandidates", "overlap", "simd", "launcher"};
	private static final String[] CATEGORIES={"CORPUS_LOWER_QUARTILE", "OTHER_CLASSIFIED", "UNKNOWN", "UNCLASSIFIED"};
	private static final String[] SUFFIXES={".all.tsv", ".frontier.tsv", ".strata.tsv", ".taxonomy.tsv", ".groups.tsv", ".gaps.tsv"};
	private static final String BREADTH_NOTE="# UNKNOWN and UNCLASSIFIED inputs count in aggregate breadth; only the rarity quantile excludes them. Whole-caller time is unmeasured.";
	private static final String HEADER="group\tfamily\tsearch\tround\tpolicy\tshipped\tinputs\t"+String.join("\t", METRICS)
		+"\trecall_pct\tprecision_pct\tfirst_gain_inputs\tlast_loss_inputs\talignment_change_vs_shipped_pct\tmodelbase_change_vs_shipped_pct"
		+"\ttp_modelbase_frontier\tbreadth_modelbase_frontier\trepeat_status\tduplicate_of\tpolicy_tuple\twhole_caller_seconds";
}
