package prok;

import java.util.ArrayList;

import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;

/** Immutable complete mapping and verifier-walk policies for one family's members.
 * Rows are contiguous by policy, with the same member order in every policy.
 * The first policy retains the shipped baseline, independently of a search's
 * moving center. Only retained names and tuples allocate during table parsing.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaJointPolicy{

	public NcrnaJointPolicy(final String name_, final Member[] members_){
		require(name_!=null && !name_.isEmpty() && members_!=null && members_.length>0,
			"A joint policy needs a name and at least one complete member tuple");
		name=name_;
		members=members_.clone();
		for(int i=0; i<members.length; i++){
			require(members[i]!=null, "Null policy member");
			for(int j=0; j<i; j++){
				require(!members[i].name.equals(members[j].name), "Duplicate policy member: "+members[i].name);
			}
		}
	}

	public int size(){return members.length;}
	public Member member(final int index){return members[index];}

	/** Parses explicit mode and all controls; fractional fixed minhits never round. */
	public static NcrnaJointPolicy[] load(final String path){
		final ArrayList<NcrnaJointPolicy> policies=new ArrayList<NcrnaJointPolicy>();
		final ArrayList<Member> members=new ArrayList<Member>();
		final ObjectIntMap<String> names=new ObjectIntMap<String>(String.class);
		final LineParser1 parser=new LineParser1('\t');
		String name=null;
		try(Input input=new Input(path)){
			final byte[] header=input.reader.nextLine();
			require(header!=null, "Missing joint policy header: "+path);
			parser.set(header);
			require(parser.terms()==COLUMNS.length, "Joint policy header must be "+HEADER);
			for(int i=0; i<COLUMNS.length; i++){
				require(parser.termEquals(COLUMNS[i], i), "Joint policy header must be "+HEADER);
			}
			for(byte[] line=input.reader.nextLine(); line!=null; line=input.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==COLUMNS.length, "Each joint policy row needs all nine columns");
				final String next=parser.parseString(0);
				if(!next.equals(name)){
					if(name!=null){append(policies, name, members);}
					require(!next.isEmpty() && names.get(next)<0, "Policy names must be nonempty and contiguous: "+next);
					names.put(next, policies.size());
					name=next;
				}
				final String member=parser.parseString(1);
				require(parser.termEquals("fixed", 2) || parser.termEquals("adaptive", 2), "Policy mode must be fixed or adaptive");
				final boolean adaptive=parser.termEquals("adaptive", 2);
				final int cap=parser.parseInt(3);
				final float floor=number(parser, 4);
				parser.setBounds(5);
				final int minHits=NcrnaShortlistTarget.fixedHits(line, parser.a(), parser.b());
				final float topFrac=number(parser, 6), qFrac=number(parser, 7);
				members.add(new Member(member, adaptive, cap, floor, minHits, topFrac, qFrac,
					Integer.parseInt(parser.parseString(8))));
			}
		}
		if(name!=null){append(policies, name, members);}
		require(!policies.isEmpty(), "A joint sweep needs a shipped baseline policy");
		return policies.toArray(new NcrnaJointPolicy[policies.size()]);
	}

	/** Configuration-only conversion preserves all finite floats, including tiny
	 * exponents and values beyond the native decimal parser's long accumulators.
	 * This allocation occurs per declared policy control, never per genomic field.
	 */
	private static float number(final LineParser1 parser, final int column){
		return Float.parseFloat(parser.parseString(column));
	}

	private static void append(final ArrayList<NcrnaJointPolicy> policies, final String name, final ArrayList<Member> members){
		final NcrnaJointPolicy point=new NcrnaJointPolicy(name, members.toArray(new Member[members.size()]));
		if(!policies.isEmpty()){
			final NcrnaJointPolicy root=policies.get(0);
			require(root.size()==point.size(), "Every policy must cover the complete baseline member roster: "+name);
			for(int i=0; i<root.size(); i++){
				require(root.member(i).name.equals(point.member(i).name), "Policy member order changed: "+name);
			}
		}
		policies.add(point);
		members.clear();
	}

	/** One member's index and patience state; omitted caller policies use packaged defaults. */
	public static final class Member{
		public Member(final String name_, final boolean adaptive_, final int topN_, final float floor_,
				final int minHits_, final float topFrac_, final float qFrac_){
			this(name_, adaptive_, topN_, floor_, minHits_, topFrac_, qFrac_, -1);
		}

		public Member(final String name_, final boolean adaptive_, final int topN_, final float floor_,
				final int minHits_, final float topFrac_, final float qFrac_, final int patience_){
			require(name_!=null && !name_.isEmpty() && topN_>0, "A member needs its registered ID and a positive cap");
			TrnaKmerIndex.validateConfiguration(1, floor_, topFrac_, qFrac_, minHits_);
			NcrnaPatience.validate(patience_);
			name=name_;
			adaptive=adaptive_;
			topN=topN_;
			floor=floor_;
			minHits=minHits_;
			topFrac=topFrac_;
			qFrac=qFrac_;
			patience=patience_;
		}

		public void apply(final NcrnaScavenger scavenger){
			scavenger.setShortlistPolicy(topN, adaptive, floor, topFrac, qFrac, minHits, patience);
		}

		/** Validate root against loaded and registered state before starting workers. */
		public void validateBaseline(final NcrnaFamily family, final String view){
			require(name.equals(family.name), "Policy roster must follow the loaded target-member order: expected "+family.name+", found "+name);
			new NcrnaShortlistTarget(name, view).validate(family, adaptive, topN, floor, minHits);
			require(patience==family.modelPatience && floor==family.adaptFloor && minHits==family.fixedMinHits
				&& topFrac==family.adaptTopFrac && qFrac==family.adaptQFrac,
				"The first policy must preserve every packaged control, including inactive values: "+name);
		}

		public final String name;
		public final boolean adaptive;
		public final int topN, minHits, patience;
		public final float floor, topFrac, qFrac;
	}

	private static void require(final boolean valid, final String why){
		if(!valid){throw new IllegalArgumentException(why);}
	}
	public final String name;
	private final Member[] members;
	private static final String[] COLUMNS={"policy", "member", "mode", "topn", "floor", "minhits", "topfrac", "qfrac", "patience"};
	public static final String HEADER="policy\tmember\tmode\ttopn\tfloor\tminhits\ttopfrac\tqfrac\tpatience";
}
