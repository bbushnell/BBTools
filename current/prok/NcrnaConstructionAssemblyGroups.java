package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

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

/** Integrates saved positive sequence identities with the original assembly groups.
 * All input rows and comparison outcomes survive. Nonmatches never join groups
 * or prove independence; a separate positive edge can still exclude an input.
 * Reads metadata only, without reopening genomes or freezing a member split.
 * @author Ganyu
 */
public final class NcrnaConstructionAssemblyGroups{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaConstructionAssemblyGroups.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String nominations=null, comparisons=null;
			int evaluation=-1, discovery=-1, pairs=-1;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("nominations")){nominations=b;}
				else if(a.equals("comparisons")){comparisons=b;}
				else if(a.equals("expectedeval")){evaluation=Parse.parseIntKMG(b);}
				else if(a.equals("expectednew")){discovery=Parse.parseIntKMG(b);}
				else if(a.equals("expectedpairs")){pairs=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown assembly-group parameter: "+arg);}
			}
			boolean assertions=false;
			assert(assertions=true);
			require(assertions && parser.in1!=null && nominations!=null && comparisons!=null && parser.out1!=null
				&& parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1,
				"Use -ea in=audit.inputs.tsv.gz nominations= comparisons= out=fresh_prefix and complete counts");
			process(parser.in1, nominations, comparisons, parser.out1, evaluation, discovery, pairs);
			pp.outstream.println("NCRNA_ASSEMBLY_GROUP_INTEGRATION_COMPLETE_NOT_FROZEN_SPLIT");
		}finally{Shared.closeStream(pp.outstream);}
	}

	public static void process(final String inputs, final String nominations, final String comparisons,
			final String output, final int evaluation, final int discovery, final int pairs){
		require(inputs!=null && nominations!=null && comparisons!=null && output!=null
			&& evaluation>0 && discovery>0 && (long)evaluation+discovery<=Integer.MAX_VALUE
			&& pairs>=0 && pairs<=1000000, "Supply complete bounded cohort and comparison counts");
		RrnaResourceIO.fresh(output, ".inputs.tsv.gz", ".pairs.tsv.gz", ".summary.tsv");
		final NcrnaConstructionAssemblyGroups audit=new NcrnaConstructionAssemblyGroups(evaluation+discovery);
		audit.load(inputs, evaluation, discovery);
		audit.join(nominations, comparisons, pairs);
		audit.write(output);
	}

	private NcrnaConstructionAssemblyGroups(final int inputs){groups=new Group[inputs];}

	/*--------------------------------------------------------------*/
	/*----------------       Original Groups        ----------------*/
	/*--------------------------------------------------------------*/

	private void load(final String path, final int expectedEvaluation, final int expectedDiscovery){
		int evaluation=0, discovery=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, INPUT_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==10 && nodes.size()<groups.length, "Original assembly audit must retain exactly ten fields and its declared population");
				final boolean isEvaluation=p.termEquals("evaluation", 0);
				require(isEvaluation || p.termEquals("discovery", 0), "Only the two original cohort namespaces are supported");
				final String id=text(row, 1), fasta=text(row, 2), canonical=text(row, 3);
				final boolean known=!p.termEquals("NA", 4);
				final int group=p.parseInt(5);
				final boolean prior=flag(6), pending=flag(7), unresolved=flag(8);
				require(!id.isEmpty() && fasta.startsWith("/") && canonical.startsWith("/")
					&& group>=0 && group<=nodes.size() && !(prior && pending),
					"Audit identities need absolute provenance paths and an earlier-or-current canonical group anchor");
				Group g=groups[group];
				if(g==null){
					require(group==nodes.size(), "The first row of an original group must be its canonical anchor");
					g=new Group(group, prior, pending, unresolved);
					groups[group]=g;
				}
				require(g.prior==prior && g.pending==pending && g.unresolvedClaim==unresolved,
					"All rows of one original group must preserve its same audit flags");
				g.containsEvaluation|=isEvaluation;
				g.known|=known;
				final Node node=new Node(row, id, fasta, isEvaluation, g, text(row, 9));
				final String key=(isEvaluation ? "evaluation\t" : "discovery\t")+id;
				require(byId.get(key)==null, "Repeated cohort-qualified input: "+key);
				byId.put(key, node);
				nodes.add(node);
				if(isEvaluation){evaluation++;}else{discovery++;}
			}
		}
		require(evaluation==expectedEvaluation && discovery==expectedDiscovery, "Original audit cohort counts differ from the declared complete population");
		for(final Node node:nodes){
			final Group g=node.group;
			require(g.containsEvaluation==g.prior, "A declared evaluation group must contain an actual evaluation input");
			final String expected=node.evaluation ? "EVALUATION" : g.prior ? KNOWN : g.pending ? PENDING
				: !g.known || g.unresolvedClaim ? UNRESOLVED : TRAIN;
			require(node.disposition.equals(expected), "Original disposition contradicts its retained group metadata: "+node.id);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Saved Pair Proof       ----------------*/
	/*--------------------------------------------------------------*/

	private void join(final String nominations, final String comparisons, final int expected){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try(Input nominated=new Input(nominations); Input compared=new Input(comparisons)){
			RrnaResourceIO.header(nominated.reader, NcrnaConstructionAssemblyCompare.HEADER, nominations);
			RrnaResourceIO.header(compared.reader, NcrnaConstructionAssemblyCompare.HEADER+"\tsequence_identity\tcriterion", comparisons);
			for(byte[] original=nominated.reader.nextLine(); original!=null; original=nominated.reader.nextLine()){
				final byte[] result=compared.reader.nextLine();
				require(pairs.size()<expected && result!=null && result.length>original.length && result[original.length]=='\t',
					"Every nominated pair needs exactly one saved result in original order");
				for(int i=0; i<original.length; i++){
					require(result[i]==original[i], "Comparison provenance differs from its original nominated row");
				}
				p.set(original);
				require(p.terms()==7 && p.termEquals(NcrnaConstructionAssemblyCompare.BASIS, 6)
					&& p.parseLong(4)>0 && p.parseLong(5)>0, "Nomination must retain positive audited dimensions and the original basis");
				final String eid=text(original, 0), did=text(original, 1);
				require(seen.add(eid+'\t'+did), "Duplicate nominated input pair would multiply saved evidence");
				final Node left=byId.get("evaluation\t"+eid), right=byId.get("discovery\t"+did);
				require(left!=null && right!=null && left.fasta.equals(text(original, 2)) && right.fasta.equals(text(original, 3))
					&& right.group.pending && !right.group.prior,
					"Each pair must bind its original evaluation/discovery IDs, literal FASTAs and pending group");
				p.set(result);
				require(p.terms()==9 && p.termEquals(NcrnaConstructionAssemblyCompare.CRITERION, 8), "Saved comparison criterion is not the reviewed complete-contig normalization");
				final boolean match=p.termEquals(NcrnaConstructionAssemblyCompare.MATCH, 7);
				require(match || p.termEquals(NcrnaConstructionAssemblyCompare.UNRESOLVED, 7), "Unknown sequence-identity outcome must not become exclusion or independence");
				right.group.nominationSeen=true;
				if(match){
					union(left.group, right.group);
					positive++;
				}else{
					root(left.group).nonmatch=true;
					root(right.group).nonmatch=true;
				}
				pairs.add(new Pair(result, left, right));
			}
			require(compared.reader.nextLine()==null && pairs.size()==expected, "Saved comparison and nomination populations must both be complete with no extra rows");
		}
		for(final Group group:groups){
			if(group==null){continue;}
			require(!group.pending || group.nominationSeen, "Every original pending group needs its saved pair evidence");
		}
	}

	private Group root(final Group group){
		Group result=group;
		while(result.parent!=result.index){result=groups[result.parent];}
		Group current=group;
		while(current.parent!=result.index){
			final Group next=groups[current.parent];
			current.parent=result.index;
			current=next;
		}
		return result;
	}

	private void union(final Group a, final Group b){
		final Group left=root(a), right=root(b);
		if(left==right){return;}
		final Group keep=left.index<right.index ? left : right;
		final Group drop=keep==left ? right : left;
		drop.parent=keep.index;
		keep.containsEvaluation|=drop.containsEvaluation;
		keep.nonmatch|=drop.nonmatch;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Complete Output       ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		long evaluation=0, known=0, added=0, pending=0, unresolved=0, train=0, groupCount=0;
		final ByteBuilder buffer=new ByteBuilder();
		try(Output inputs=new Output(prefix+".inputs.tsv.gz"); Output links=new Output(prefix+".pairs.tsv.gz");
				Output summary=new Output(prefix+".summary.tsv")){
			inputs.writer.println(INPUT_HEADER+"\tresolved_group_id\tresolved_evaluation_group\tgroup_has_unresolved_pair\tresolved_assembly_disposition");
			links.writer.println(NcrnaConstructionAssemblyCompare.HEADER+"\tsequence_identity\tcriterion\tresolved_evaluation_group_id\tresolved_discovery_group_id");
			for(final Node node:nodes){
				final Group group=root(node.group);
				final String disposition;
				if(node.evaluation){
					disposition="EVALUATION";
					evaluation++;
				}else if(node.group.prior){
					disposition=KNOWN;
					known++;
				}else if(group.containsEvaluation){
					disposition=VERIFIED;
					added++;
				}else if(node.group.pending){
					disposition=SEQUENCE_UNRESOLVED;
					pending++;
				}else if(node.disposition.equals(UNRESOLVED)){
					disposition=UNRESOLVED;
					unresolved++;
				}else{
					disposition=TRAIN;
					train++;
				}
				inputs.writer.print(buffer.clear().append(node.row).tab().append(group.index).tab()
					.append(group.containsEvaluation).tab().append(group.nonmatch).tab().append(disposition).nl());
			}
			for(final Pair pair:pairs){
				links.writer.print(buffer.clear().append(pair.row).tab().append(root(pair.left.group).index)
					.tab().append(root(pair.right.group).index).nl());
			}
			for(final Group group:groups){if(group!=null && root(group)==group){groupCount++;}}
			require(evaluation+known+added+pending+unresolved+train==nodes.size(), "Integrated disjoint dispositions must conserve every original input");
			summary.writer.println("evaluation_inputs\tknown_excluded_discovery\tnew_sequence_excluded_discovery\tsequence_unresolved_discovery\tother_unresolved_discovery\ttraining_candidate_discovery\texact_match_pairs\tunresolved_pairs\tresolved_groups\tstatus");
			summary.writer.print(buffer.clear().append(evaluation).tab().append(known).tab().append(added).tab()
				.append(pending).tab().append(unresolved).tab().append(train).tab().append(positive).tab()
				.append(pairs.size()-positive).tab().append(groupCount).append("\tASSEMBLY_GROUPS_NOT_FROZEN_MEMBER_SPLIT\n"));
		}
	}

	private boolean flag(final int column){
		require(p.termEquals("true", column) || p.termEquals("false", column), "Audit flags must be literal true or false");
		return p.termEquals("true", column);
	}
	private String text(final byte[] row, final int column){
		p.setBounds(column);
		final String value=new String(row, p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\ufffd')<0, "Invalid UTF-8 must not change retained provenance");
		return value;
	}
	private static final class Group{
		Group(final int index_, final boolean prior_, final boolean pending_, final boolean unresolved_){
			index=index_;
			parent=index;
			prior=prior_;
			pending=pending_;
			unresolvedClaim=unresolved_;
		}
		final int index;
		final boolean prior, pending, unresolvedClaim;
		int parent;
		boolean containsEvaluation, known, nominationSeen, nonmatch;
	}
	private static final class Node{
		Node(final byte[] row_, final String id_, final String fasta_, final boolean evaluation_, final Group group_, final String disposition_){
			row=row_;
			id=id_;
			fasta=fasta_;
			evaluation=evaluation_;
			group=group_;
			disposition=disposition_;
		}
		final byte[] row;
		final String id, fasta, disposition;
		final boolean evaluation;
		final Group group;
	}
	private static final class Pair{
		Pair(final byte[] row_, final Node left_, final Node right_){
			row=row_;
			left=left_;
			right=right_;
		}
		final byte[] row;
		final Node left, right;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final Group[] groups;
	private final ArrayList<Node> nodes=new ArrayList<Node>();
	private final ArrayList<Pair> pairs=new ArrayList<Pair>();
	private final ObjectMap<String, Node> byId=new ObjectMap<String, Node>(String.class, Node.class);
	private final LineParser1 p=new LineParser1('\t');
	private long positive;
	public static final String INPUT_HEADER="cohort\tinput_id\tfasta\tcanonical_fasta\tassembly_record\tgroup_id\tevaluation_group\tsequence_comparison_pending\tunresolved_claim\tdisposition";
	private static final String KNOWN="EXCLUDE_KNOWN_EVALUATION_ASSEMBLY";
	private static final String VERIFIED="EXCLUDE_VERIFIED_EVALUATION_ASSEMBLY_ALIAS";
	private static final String PENDING="ASSEMBLY_SEQUENCE_CHECK_PENDING";
	private static final String SEQUENCE_UNRESOLVED="ASSEMBLY_SEQUENCE_IDENTITY_UNRESOLVED";
	private static final String UNRESOLVED="ASSEMBLY_IDENTITY_UNRESOLVED";
	private static final String TRAIN="TRAIN_CANDIDATE_ASSEMBLY_IDENTITY";
}
