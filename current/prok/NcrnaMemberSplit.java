package prok;

import java.util.ArrayList;
import java.util.Locale;

import map.IntHashMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Rebinds saved FMN admission rows to the common assembly split.
 * Admission, native witnesses, original cohort and sequence-validation status
 * remain intact. Metadata candidate counts are not sequence-deduplicated weights.
 * @author Ganyu
 */
public final class NcrnaMemberSplit{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaMemberSplit.class, false);
		try{
			final Config c=new Config();
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String a=pair[0].toLowerCase(Locale.ROOT);
				final String b=pair.length>1 ? pair[1] : null;
				if(a.equals("split")){c.split=b;}
				else if(a.equals("groups")){c.groups=b;}
				else if(a.equals("expectedinputs")){c.inputs=Parse.parseIntKMG(b);}
				else if(a.equals("expectedpool")){c.pool=Parse.parseIntKMG(b);}
				else if(a.equals("expectedold")){c.old=Parse.parseIntKMG(b);}
				else if(a.equals("expectedevents")){c.events=Parse.parseIntKMG(b);}
				else if(a.equals("expectedcores")){c.cores=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown member-split option: "+arg);}
			}
			c.members=parser.in1;
			c.output=parser.out1;
			require(parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1, "Reassign complete saved ledgers without sampling or overwrite");
			run(c);
			pp.outstream.println("NCRNA_MEMBER_SPLIT_COMPLETE_GROUP_SEQUENCE_DEDUP_REQUIRED");
		}finally{Shared.closeStream(pp.outstream);}
	}
	public static void run(final Config c){
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Ledger rebinding requires native assertions");
		require(c!=null && c.members!=null && c.split!=null && c.groups!=null && c.output!=null
			&& c.inputs>0 && c.pool>0 && c.old>=0 && c.events>=0 && c.cores>=0, "Declare all saved inputs and expected row populations");
		RrnaResourceIO.fresh(c.output, ".old.tsv.gz", ".events.tsv.gz", ".extract.tsv.gz", ".labels.tsv", ".summary.tsv");
		final NcrnaMemberSplit x=new NcrnaMemberSplit();
		x.nodes(c.groups, c.inputs);
		x.split(c.split, c.pool);
		x.old(c.members+".old.tsv.gz", c.output, c.old);
		x.events(c.members+".events.tsv.gz", c.output, c.events);
		x.cores(c.members+".extract.tsv.gz", c.output, c.cores);
		x.links();
		x.summary(c);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	private void nodes(final String path, final int expected){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.GROUP_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 14);
				final Node n=new Node(text(0), p.parseInt(1), text(2), p.parseInt(10));
				n.unresolved=flag(12);
				require(n.id>0 && n.group>=0 && n.path.startsWith("/") && (n.cohort.equals("evaluation") || n.cohort.equals("discovery"))
					&& byInput.get(n.key())==null && byPath.get(n.path)==null, "Original resolved input identities must be complete and unique");
				byInput.put(n.key(), n);
				byPath.put(n.path, n);
				nodes.add(n);
				resolvedGroups.put(n.group, 1);
			}
		}
		require(nodes.size()==expected, "All original evaluation/discovery inputs must survive label rebinding");
	}
	private void split(final String prefix, final int expected){
		try(Input in=new Input(prefix+".groups.tsv")){
			RrnaResourceIO.header(in.reader, "assembly_group\tstratum\tpool_inputs\tassignment", prefix);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 4);
				final int group=p.parseInt(0);
				final int label=label(text(3));
				require(resolvedGroups.contains(group) && !groupIndex.contains(group) && label<OUTSIDE && p.parseInt(2)>0,
					"Common split groups must be unique original groups with positive pool membership");
				groupIndex.put(group, groups.size());
				groups.add(new Group(group, text(1), label, p.parseInt(2)));
			}
		}
		int count=0;
		try(Input in=new Input(prefix+".inputs.tsv.gz")){
			RrnaResourceIO.header(in.reader, "assembly\tfasta\tassembly_group\tstratum\tassignment", prefix);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 5);
				final Node n=byPath.get(text(1));
				final Group g=group(p.parseInt(2));
				require(n!=null && !n.inPool && n.group==g.id && p.termEquals(g.stratum, 3) && label(text(4))==g.label,
					"Every pool input must bind its original path/group and common stratum/side");
				n.inPool=true;
				g.observed++;
				count++;
			}
		}
		require(count==expected, "The complete common PGAP pool must be present");
		for(final Group g:groups){require(g.observed==g.expected, "Group membership must conserve each declared common-split quota unit");}
		final ObjectSet<String> pairPaths=new ObjectSet<String>(String.class);
		int pairAnchor=-1, pairPartner=-1;
		try(Input in=new Input(prefix+".coassigned.tsv")){
			RrnaResourceIO.header(in.reader, "pool_group\tpartner_group\tinput_group\tfasta\tin_PGAP_pool\tassignment\treason", prefix);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 7);
				final Group anchor=group(p.parseInt(0));
				final int partner=p.parseInt(1);
				if(pairAnchor<0){pairAnchor=anchor.id; pairPartner=partner;}
				require(anchor.id==pairAnchor && partner==pairPartner && anchor.id!=partner,
					"The common split contains one reviewed unresolved pair, not inferred extra exceptions");
				final Node n=byPath.get(text(3));
				require(n!=null && n.group==p.parseInt(2) && (n.group==anchor.id || n.group==partner)
					&& flag(4)==n.inPool && label(text(5))==anchor.label
					&& p.termEquals("unresolved, co-assigned to prevent leakage", 6) && pairPaths.add(n.path),
					"Explicit pair provenance must retain each original input/group, pool membership and same side exactly once");
				final Group current=group(n.group);
				require(current.label==anchor.label && current.stratum.equals(anchor.stratum), "Co-assignment cannot override a common pool assignment or stratum");
				current.reason="unresolved, co-assigned to prevent leakage";
			}
		}
		for(final Node n:nodes){
			if(n.group==pairAnchor || n.group==pairPartner){require(pairPaths.contains(n.path), "Unresolved pair provenance must retain every original alias path");}
			require(!n.inPool || !n.unresolved || !reason(n.group).equals("NA"), "An unresolved pool group needs its explicit accepted pair provenance");
			n.label=assignment(n.group);
		}
	}
	private void old(final String path, final String output, final int expected){
		try(Input in=new Input(path); Output out=new Output(output+".old.tsv.gz")){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.OLD_OUTPUT_HEADER, path);
			out.writer.println(OLD_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 39);
				final Old old=new Old(text(5), p.parseInt(34), flag(37), flag(38));
				require(resolvedGroups.contains(old.group) && oldById.get(old.id)==null && (!old.intact || old.seed)
					&& (!flag(36) || !old.intact && !old.seed), "Saved old-member eligibility and contradiction flags must remain coherent");
				oldById.put(old.id, old);
				olds.add(old);
				old.coreIds=text(35);
				final int label=assignment(old.group);
				final boolean train=label!=HELD;
				out.writer.print(b.clear().append(row).tab().append(LABELS[label]).tab()
					.append(train && old.intact).tab().append(train && old.seed).tab().append(reason(old.group)).nl());
				final long[] sums=counts[label];
				sums[0]++;
				if(old.intact){sums[1]++;}
				if(old.seed){sums[2]++;}
				if(train && old.intact){sums[3]++;}
				if(train && old.seed){sums[4]++;}
			}
		}
		require(olds.size()==expected, "Every existing old-member row must survive the new split");
	}
	private void events(final String path, final String output, final int expected){
		final String header=NcrnaConstructionMembers.FMN_HEADER+"\tadmission\tassembly_split\tassembly_group\textraction_member";
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try(Input in=new Input(path); Output out=new Output(output+".events.tsv.gz")){
			RrnaResourceIO.header(in.reader, header, path);
			out.writer.println(EVENT_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 36);
				final Node n=node(text(0), p.parseInt(3));
				require(n.group==p.parseInt(34) && p.termEquals("RF00050", 8) && ids.add(text(0)+'\t'+text(2)),
					"Native events retain their original input/group/family and unique identity");
				final String core=text(35);
				if(!core.equals("NA")){eventLinks.add(new Link(core, n.group));}
				out.writer.print(b.clear().append(row).tab().append(LABELS[n.label]).tab().append(reason(n.group)).nl());
				counts[n.label][5]++;
			}
		}
		require(ids.size()==expected, "Every saved native event must survive, including nonadmitted events");
	}
	private void cores(final String path, final String output, final int expected){
		try(Input in=new Input(path); Output out=new Output(output+".extract.tsv.gz")){
			RrnaResourceIO.header(in.reader, NcrnaConstructionMembers.EXTRACT_HEADER, path);
			out.writer.println(EXTRACT_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				parse(row, 21);
				final Node n=node(text(1), p.parseInt(2));
				require(n.group==p.parseInt(4) && p.termEquals(n.path, 3) && p.parseInt(11)>0 && coreGroups.get(text(0))==null,
					"Each retained extraction core must keep its original source path/group and unique identifier");
				coreGroups.put(text(0), n);
				final boolean intact=flag(14), seed=flag(15);
				require(!intact || seed, "An intact admitted core remains eligible for seed construction");
				boolean oldIntact=false, oldSeed=false;
				final String ids=text(18);
				if(!ids.equals("NA")){
					int start=0;
					while(start<ids.length()){
						int end=ids.indexOf(';', start);
						if(end<0){end=ids.length();}
						final Old old=oldById.get(ids.substring(start, end));
						require(old!=null && old.group==n.group, "A core's old-member associations must preserve resolved assembly identity");
						oldIntact|=old.intact;
						oldSeed|=old.seed;
						start=end+1;
					}
				}
				final boolean trainIntact=intact && (n.label==TRAIN || n.label==CRITICAL || n.label==OUTSIDE && oldIntact);
				final boolean trainSeed=seed && (n.label==TRAIN || n.label==CRITICAL || n.label==OUTSIDE && oldSeed);
				require(!trainIntact || trainSeed, "Methodology intact eligibility cannot exceed seed eligibility");
				out.writer.print(b.clear().append(row).tab().append(LABELS[n.label]).tab().append(trainIntact)
					.tab().append(trainSeed).tab().append(reason(n.group)).append("\tGROUP_SEQUENCE_DEDUP_PENDING\n"));
				final long[] sums=counts[n.label];
				sums[6]++;
				if(intact){sums[7]++;}
				if(seed){sums[8]++;}
				if(trainIntact){sums[9]++;}
				if(trainSeed){sums[10]++;}
			}
		}
		require(coreGroups.size()==expected, "Every saved candidate core must remain in the ledger regardless of its new side");
	}
	private void links(){
		for(final Link link:eventLinks){require(coreGroups.get(link.id)!=null && coreGroups.get(link.id).group==link.group, "Every admitted event must still point to its retained core/group");}
		for(final Old old:olds){
			if(old.coreIds.equals("NA")){continue;}
			for(final String id:old.coreIds.split(";")){
				final Node node=coreGroups.get(id);
				require(node!=null && node.group==old.group, "Every old-member extraction link must retain its original resolved group");
			}
		}
	}
	private void summary(final Config c){
		long old=0, events=0, cores=0;
		try(Output out=new Output(c.output+".labels.tsv"); Output summary=new Output(c.output+".summary.tsv")){
			out.writer.println("assignment\told_rows\told_intact_eligible\told_seed_eligible\told_methodology_intact\told_methodology_seed\tnative_events\tcandidate_core_rows\tintact_candidate_rows\tseed_candidate_rows\tmethodology_intact_core_rows\tmethodology_seed_core_rows");
			for(int label=0; label<LABELS.length; label++){
				b.clear().append(LABELS[label]);
				for(final long count:counts[label]){b.tab().append(count);}
				out.writer.print(b.nl());
				old+=counts[label][0];
				events+=counts[label][5];
				cores+=counts[label][6];
			}
			require(old==c.old && events==c.events && cores==c.cores, "Disjoint assignment labels must conserve every existing ledger row");
			summary.writer.println("old_members\tnative_events\tcandidate_core_rows\tcommon_pool_inputs\tresolved_inputs\tstatus");
			summary.writer.print(b.clear().append(old).tab().append(events).tab().append(cores).tab().append(c.pool).tab().append(c.inputs)
				.append("\tMETADATA_REBOUND_GROUP_SEQUENCE_DEDUP_REQUIRED_BEFORE_LEARNING\n"));
		}
	}
	private Group group(final int id){final int i=groupIndex.get(id); require(i>=0, "Common split input refers to an absent group"); return groups.get(i);}
	private int assignment(final int id){final int i=groupIndex.get(id); return i<0 ? OUTSIDE : groups.get(i).label;}
	private String reason(final int id){final int i=groupIndex.get(id); return i<0 ? "NA" : groups.get(i).reason;}
	private Node node(final String cohort, final int id){final Node n=byInput.get(cohort+'\t'+id); require(n!=null, "Saved member/event input must bind the original resolved catalog"); return n;}
	private static int label(final String value){for(int i=0; i<LABELS.length; i++){if(LABELS[i].equals(value)){return i;}} throw new IllegalArgumentException("Unknown common assignment: "+value);}
	private void parse(final byte[] row, final int fields){p.set(row); require(p.terms()==fields, "Saved member-split schema changed");}
	private String text(final int field){return p.parseString(field);}
	private boolean flag(final int field){require(p.termEquals("true", field) || p.termEquals("false", field), "Eligibility must be an explicit boolean"); return p.termEquals("true", field);}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String members, split, groups, output;
		public int inputs=41082, pool=39979, old=26896, events=27666, cores=27631;
	}
	private static final class Node{
		Node(final String c, final int i, final String p, final int g){
			cohort=c;
			id=i;
			path=p;
			group=g;
		}
		String key(){return cohort+'\t'+id;}
		final String cohort, path;
		final int id, group;
		boolean inPool, unresolved;
		int label;
	}
	private static final class Group{
		Group(final int i, final String s, final int l, final int n){
			id=i;
			stratum=s;
			label=l;
			expected=n;
		}
		final int id, label, expected;
		final String stratum;
		int observed;
		String reason="NA";
	}
	private static final class Old{
		Old(final String i, final int g, final boolean a, final boolean s){
			id=i;
			group=g;
			intact=a;
			seed=s;
		}
		final String id;
		final int group;
		final boolean intact, seed;
		String coreIds;
	}
	private static final class Link{
		Link(final String i, final int g){
			id=i;
			group=g;
		}
		final String id;
		final int group;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final LineParser1 p=new LineParser1('\t');
	private final ByteBuilder b=new ByteBuilder();
	private final ObjectMap<String, Node> byInput=new ObjectMap<String, Node>(String.class, Node.class), byPath=new ObjectMap<String, Node>(String.class, Node.class), coreGroups=new ObjectMap<String, Node>(String.class, Node.class);
	private final ObjectMap<String, Old> oldById=new ObjectMap<String, Old>(String.class, Old.class);
	private final IntHashMap resolvedGroups=new IntHashMap(), groupIndex=new IntHashMap();
	private final ArrayList<Node> nodes=new ArrayList<Node>();
	private final ArrayList<Group> groups=new ArrayList<Group>();
	private final ArrayList<Old> olds=new ArrayList<Old>();
	private final ArrayList<Link> eventLinks=new ArrayList<Link>();
	private final long[][] counts=new long[4][11];
	private static final int TRAIN=0, CRITICAL=1, HELD=2, OUTSIDE=3;
	private static final String[] LABELS={"TRAIN", "COVERAGE_CRITICAL_TRAIN", "HELD_OUT", "OUTSIDE_POOL"};
	public static final String OLD_HEADER=NcrnaConstructionMembers.OLD_OUTPUT_HEADER+"\tmethodology_assignment\tmethodology_intact_eligible\tmethodology_seed_eligible\tcoassignment_reason";
	public static final String EVENT_HEADER=NcrnaConstructionMembers.FMN_HEADER+"\tadmission\tassembly_split\tassembly_group\textraction_member\tmethodology_assignment\tcoassignment_reason";
	public static final String EXTRACT_HEADER=NcrnaConstructionMembers.EXTRACT_HEADER+"\tmethodology_assignment\tmethodology_train_intact\tmethodology_train_seed\tcoassignment_reason\tweight_status";
}
