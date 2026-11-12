package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

import map.IntObjectMap;
import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaReferenceCalls.Rule;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import structures.ByteBuilder;

/** Reads retained member evidence without reclassifying native scores. The grading
 * constructor retains the cobalamin-only contract; inventory supplies an explicit
 * family/member/model map and never infers membership from caller prefixes.
 * One thread owns the readers and selected-parent index. Both sidecars are streamed
 * once; only selected inputs and mapped models are retained. Source identity must
 * be authenticated by the invoking project against its registered manifest.
 * Binding checks here establish consistency, not biological truth or authorization.
 * Protected prior assignments override overlapping current-model observations.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaSubtypeEvidence{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	NcrnaSubtypeEvidence(final String mapping, final String witnesses, final String events,
		final Rule[] rules, final List<String> inputs){
		require(inputs!=null && rules!=null, "Subtype evidence needs the complete selected input population and query rules");
		witnessPath=witnesses;
		for(final String input:inputs){selectedInputs.put(input, input);}
		loadTypes(mapping, rules);
		final IntObjectMap<Event> byEvent=loadEvents(events);
		loadWitnesses(witnesses, byEvent);
	}

	/** Metadata inventory has no caller queries; model membership is explicit. */
	NcrnaSubtypeEvidence(final String mapping, final String witnesses, final String events,
		final List<String> inputs){
		require(inputs!=null && !inputs.isEmpty(), "Inventory support needs the complete original input population");
		witnessPath=witnesses;
		for(final String input:inputs){
			require(!selectedInputs.contains(input), "Repeated inventory input would duplicate support");
			selectedInputs.put(input, input);
		}
		try(Input in=new Input(mapping)){
			RrnaResourceIO.header(in.reader, MEMBER_HEADER, mapping);
			for(byte[] line; (line=in.reader.nextLine())!=null;){
				parser.set(line);
				require(parser.terms()==3, "Member mapping requires family, member and reference_model");
				final Type type=new Type(text(1), null, text(2), text(0));
				require(!type.family.isEmpty() && !type.name.isEmpty() && !type.model.isEmpty()
					&& !type.name.equals(".") && type.name.indexOf(',')<0,
					"Member fields must be nonempty; names cannot collide with support-list markers");
				require(types.size()<30, "Member support uses at most thirty positive bit positions");
				for(final Type prior:types){
					require(!prior.model.equals(type.model)
						&& !(prior.family.equals(type.family) && prior.name.equals(type.name)),
						"Reference models and family/member pairs must be unique");
				}
				types.add(type);
			}
		}
		require(!types.isEmpty(), "Member mapping cannot be empty");
		loadWitnesses(witnesses, loadEvents(events));
	}

	/** Mapping is explicit and complete for included cobalamin query prefixes. */
	private void loadTypes(final String path, final Rule[] rules){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, MAP_HEADER, path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==3, "Subtype mapping requires subtype, model_prefix and reference_model");
				final Type type=new Type(text(0), text(1), text(2));
				require(!type.name.isEmpty() && !type.prefix.isEmpty() && !type.model.isEmpty(), "Subtype mapping fields cannot be empty");
				require(!type.name.equals(".") && type.name.indexOf(',')<0,
					"Subtype names cannot collide with the missing marker or support-list separator");
				require(types.size()<30, "Subtype support uses at most thirty positive bit positions");
				for(final Type prior:types){
					require(!prior.name.equals(type.name) && !prior.prefix.equals(type.prefix) && !prior.model.equals(type.model),
						"Subtype names, query prefixes and reference models must each be unique");
				}
				int found=-1;
				for(int i=0; i<rules.length; i++){
					if(rules[i].modelPrefix.equals(type.prefix)){
						require(FAMILY.equals(rules[i].family) && rules[i].include, "Subtype query prefix must map to included cobalamin");
						found=i;
					}
				}
				require(found>=0, "Subtype prefix is absent from the explicit query family map: "+type.prefix);
				type.rule=found;
				types.add(type);
			}
		}
		require(!types.isEmpty(), "Subtype mapping cannot be empty");
		queryTypes=new int[rules.length];
		java.util.Arrays.fill(queryTypes, -1);
		for(int i=0; i<types.size(); i++){queryTypes[types.get(i).rule]=i;}
		for(int i=0; i<rules.length; i++){
			if(FAMILY.equals(rules[i].family) && rules[i].include){
				require(queryTypes[i]>=0, "Every included cobalamin query prefix requires subtype mapping: "+rules[i].modelPrefix);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Retained Evidence      ----------------*/
	/*--------------------------------------------------------------*/

	/** Keeps native event identities and already-decided tiers, never recomputes GA. */
	private IntObjectMap<Event> loadEvents(final String path){
		final IntObjectMap<Event> events=new IntObjectMap<Event>(Event.class);
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, EVENT_HEADER, path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==30, "Subtype event evidence requires the thirty-field tier provenance schema");
				if(!mappedFamily(7)){continue;}
				final String input=selectedInput(1);
				if(input==null){continue;}
				final int type=model(6);
				require(type>=0 && parser.termEquals(types.get(type).family, 7),
					"Selected family event has an unmapped or misbound reference model");
				int tier=-1;
				for(int i=0; i<TIERS.length; i++){if(parser.termEquals(TIERS[i], 22)){tier=i;}}
				require(tier>=0, "Unknown retained event tier cannot confer subtype support");
				// NcrnaReferenceTiers.groups excludes INVALID_INTERVAL from every witness group.
				if(tier==INVALID){continue;}
				final int id=positiveInt(0), start=positiveInt(3), stop=positiveInt(4);
				require(stop>=start, "Selected subtype event interval must be ordered");
				final char strand=parser.termEquals("+", 5) ? '+' : parser.termEquals("-", 5) ? '-' : '.';
				require(strand!='.', "Selected subtype event requires an explicit strand");
				final long nativeLine=parser.parseLong(27);
				require(nativeLine>0 && !events.contains(id), "Subtype event IDs must be unique and native lines positive");
				events.put(id, new Event(input, decoded(2), start, stop, strand, type, tier,
					(parser.parseInt(25)&BRIDGE_FLAG)!=0, text(26), nativeLine));
			}
		}
		return events;
	}

	/** Binds original witness records to retained events before aggregating support. */
	private void loadWitnesses(final String path, final IntObjectMap<Event> events){
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, WITNESS_HEADER, path);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				parser.set(line);
				require(parser.terms()==11, "Annotation witnesses require eleven fields");
				final int type=model(6);
				if(type<0){continue;}
				final String input=selectedInput(1);
				if(input==null){continue;}
				final Support support=parent(input);
				final String id=support.id;
				require(support.input.equals(input), "A parent witness ID cannot belong to two inputs", id);
				final int start=positiveInt(7), stop=positiveInt(8);
				final long priorLine=parser.parseLong(9);
				require(stop>=start, "Subtype witness interval must be ordered");
				if(parser.termEquals("PROTECTED_PRIOR", 2)){
					require(support.priorType<0 && priorLine>0 && parser.termEquals(".", 3)
						&& parser.termEquals(".", 4) && parser.termEquals(".", 5)
						&& parser.termEquals("ORIGINAL_ACCEPTED_INTERVAL_AND_AUTHORITY", 10),
						"A protected parent requires exactly one original prior authority witness");
					support.priorType=type;
					support.priorLine=priorLine;
					support.priorStart=start;
					support.priorStop=stop;
				}else{
					require(parser.termEquals("CM_HIT", 2), "Unknown annotation evidence kind");
					final Event event=events.get(positiveInt(3));
					require(event!=null && event.input.equals(input) && event.type==type
						&& event.start==start && event.stop==stop && event.line==parser.parseLong(5)
						&& parser.termEquals(event.table, 4), "Witness event/input/model/interval/native binding mismatch", id);
					int association=-1;
					for(int i=0; i<ASSOCIATIONS.length; i++){if(parser.termEquals(ASSOCIATIONS[i], 10)){association=i;}}
					require(association>=0, "Unknown subtype witness association", id);
					require(association==PRIOR_OVERLAP ? priorLine>0 : priorLine==-1,
						"Prior-line association disagrees with witness role", id);
					if(association==GROUP){
						require(!event.bridge && (event.tier==CONFIDENT || event.tier==MARGINAL),
							"Grouping subtype support requires a non-bridge confident or marginal event");
					}else if(association==BRIDGE){
						require(event.bridge, "Tandem-bridge witness must retain its saved bridge flag");
					}else if(association==WEAK){
						require(!event.bridge && event.tier!=CONFIDENT && event.tier!=MARGINAL,
							"Weak witness cannot replace qualifying or bridge provenance");
					}
					for(final Link previous:support.links){require(previous.event!=event, "Repeated event witness for one parent", id);}
					support.links.add(new Link(event, association, priorLine));
				}
			}
		}
	}

	/** Consecutive witnesses share a parent in the tier writer. Reuse its decoded
	 * identity without allocating a String per event; the map also accepts later repeats. */
	private Support parent(final String input){
		parser.setBounds(0);
		final int from=parser.a(), length=parser.b()-from;
		boolean same=lastParent!=null && parentToken.length()==length;
		for(int i=0; same && i<length; i++){same=parentToken.array[i]==parser.line()[from+i];}
		if(same){return lastParent;}
		final String id=decoder.decode(parser.line(), from, from+length);
		Support support=parents.get(id);
		if(support==null){
			support=new Support(id, input);
			parents.put(id, support);
			order.add(support);
		}
		parentToken.clear().append(parser.line(), from, length);
		lastParent=support;
		return support;
	}

	/** Checks a selected authoritative parent, then derives distinct subtype bits.
	 * No interval is changed and weak/bridge observations never create masks. */
	Support support(final Row row){
		require(hasFamily(row.locus.family), "Support projection requires a declared reference family");
		final Support support=parents.get(row.id);
		require(support!=null && support.input.equals(row.input) && !support.used,
			"Missing, repeated or misbound subtype parent evidence", row.id);
		require(witnessPath.equals(row.evidence), "Subtype parent evidence path differs from the supplied witness table", row.id);
		for(final Link link:support.links){
			final Event event=link.event;
			require(types.get(event.type).family.equals(row.locus.family), "Witness member belongs to another reference family", row.id);
			require(event.target.equals(row.locus.seqid) && event.strand==row.locus.strand
				&& event.start<=row.locus.stop && row.locus.start<=event.stop, "Subtype witness does not overlap its parent on the same contig/strand", row.id);
			if(support.priorType>=0){
				require(link.association==PRIOR_OVERLAP && link.priorLine==support.priorLine,
					"Protected prior cannot acquire regrouped subtype authority", row.id);
			}else{
				require(link.association!=PRIOR_OVERLAP, "Prior-associated event lacks protected authority", row.id);
				if(link.association==GROUP){
					if(event.tier==CONFIDENT){support.scored|=1<<event.type;}
					else{support.marginal|=1<<event.type;}
				}
			}
		}
		if(support.priorType>=0){
			require(!row.mask && types.get(support.priorType).family.equals(row.locus.family)
				&& types.get(support.priorType).model.equals(row.model)
				&& row.locus.start==support.priorStart && row.locus.stop==support.priorStop,
				"Protected subtype assignment and interval must retain their original authority", row.id);
			support.scored=1<<support.priorType;
		}
		require(row.mask || support.scored!=0, "Scored parent lacks justified scored subtype support", row.id);
		require(!row.mask || support.scored==0, "Masked parent cannot carry confident scored subtype support", row.id);
		int representative=-1;
		for(int i=0; i<types.size(); i++){if(types.get(i).model.equals(row.model)){representative=i;}}
		require(representative>=0 && ((support.scored|support.marginal)&(1<<representative))!=0,
			"Parent representative lacks its original qualifying subtype witness", row.id);
		support.marginal&=~support.scored;
		support.used=true;
		return support;
	}

	/** Detects witness parents omitted from the selected reference population. */
	void finish(){
		for(final Support support:order){require(support.used, "Selected subtype witness has no reference parent", support.id);}
	}
	int queryType(final int rule){return queryTypes[rule];}
	int size(){return types.size();}
	String name(final int index){return types.get(index).name;}
	String family(final int index){return types.get(index).family;}
	String referenceModel(final int index){return types.get(index).model;}
	boolean hasFamily(final String family){
		for(final Type type:types){if(type.family.equals(family)){return true;}}
		return false;
	}
	private boolean mappedFamily(final int column){
		for(final Type type:types){if(parser.termEquals(type.family, column)){return true;}}
		return false;
	}

	/** Shared with the cobalamin grader: masks survive in every member view. */
	static boolean includedInView(final Row row, final Support item, final int type){
		assert(row!=null && item!=null && type>=0 && type<30) : "View membership needs a bound parent and valid member bit";
		return row.mask || ((item.scored|item.marginal)&(1<<type))!=0;
	}
	/** Exclusive support masks a multi-supported parent before matching, never after. */
	static boolean maskedInView(final Row row, final Support item, final int type, final boolean exclusive){
		assert(includedInView(row, item, type)) : "Only included reference parents have a view disposition";
		return row.mask || (item.scored&(1<<type))==0 || exclusive && item.multi();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Parsing           ----------------*/
	/*--------------------------------------------------------------*/

	private int model(final int column){
		for(int i=0; i<types.size(); i++){if(parser.termEquals(types.get(i).model, column)){return i;}}
		return -1;
	}
	/** Input IDs commonly repeat; allocate on input transitions, not every hit. */
	private String selectedInput(final int column){
		if(lastInput==null || !parser.termEquals(lastInput, column)){
			lastInput=text(column);
			lastSelected=selectedInputs.get(lastInput);
		}
		return lastSelected;
	}
	private int positiveInt(final int column){
		final long value=parser.parseLong(column);
		require(value>0 && value<=Integer.MAX_VALUE, "Event IDs and coordinates must fit positive integer indexing");
		return (int)value;
	}
	private String decoded(final int column){
		parser.setBounds(column);
		return decoder.decode(parser.line(), parser.a(), parser.b());
	}
	private String text(final int column){
		parser.setBounds(column);
		final String value=new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
		require(value.indexOf('\uFFFD')<0, "Invalid UTF-8 would change subtype evidence identity");
		return value;
	}
	private static void require(final boolean valid, final String reason){
		if(!valid){throw new IllegalArgumentException(reason);}
	}
	/** Contextual messages allocate only on failure, not once per witness. */
	private static void require(final boolean valid, final String reason, final String parent){
		if(!valid){throw new IllegalArgumentException(reason+": "+parent);}
	}

	/*--------------------------------------------------------------*/
	/*----------------             Types            ----------------*/
	/*--------------------------------------------------------------*/

	static final class Support{
		Support(final String id_, final String input_){
			id=id_;
			input=input_;
		}
		boolean multi(){return Integer.bitCount(scored)>1;}
		int bridgeWitnesses(){
			int count=0;
			for(final Link link:links){if(link.association==BRIDGE){count++;}}
			return count;
		}
		final String id, input;
		final ArrayList<Link> links=new ArrayList<Link>();
		int scored, marginal, priorType=-1, priorStart, priorStop;
		long priorLine;
		boolean used;
	}
	private static final class Type{
		Type(final String name_, final String prefix_, final String model_){
			this(name_, prefix_, model_, FAMILY);
		}
		Type(final String name_, final String prefix_, final String model_, final String family_){
			name=name_;
			prefix=prefix_;
			model=model_;
			family=family_;
		}
		final String name, prefix, model, family;
		int rule;
	}
	private static final class Event{
		Event(final String input_, final String target_, final int start_, final int stop_, final char strand_,
			final int type_, final int tier_, final boolean bridge_, final String table_, final long line_){
			input=input_;
			target=target_;
			start=start_;
			stop=stop_;
			strand=strand_;
			type=type_;
			tier=tier_;
			bridge=bridge_;
			table=table_;
			line=line_;
		}
		final String input, target, table;
		final int start, stop, type, tier;
		final char strand;
		final long line;
		final boolean bridge;
	}
	private static final class Link{
		Link(final Event event_, final int association_, final long priorLine_){
			event=event_;
			association=association_;
			priorLine=priorLine_;
		}
		final Event event;
		final int association;
		final long priorLine;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String witnessPath;
	private final LineParser1 parser=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final ObjectMap<String, String> selectedInputs=new ObjectMap<String, String>(String.class, String.class);
	private final ObjectMap<String, Support> parents=new ObjectMap<String, Support>(String.class, Support.class);
	private final ArrayList<Support> order=new ArrayList<Support>();
	private final ArrayList<Type> types=new ArrayList<Type>();
	private String lastInput, lastSelected;
	private Support lastParent;
	private final ByteBuilder parentToken=new ByteBuilder();
	private int[] queryTypes;
	static final String FAMILY="cobalamin";
	static final String MAP_HEADER="subtype\tmodel_prefix\treference_model";
	static final String MEMBER_HEADER="family\tmember\treference_model";
	static final String WITNESS_HEADER="annotation_id\tinput_id\tevidence_kind\tevent_id\tnative_table\tnative_line\tmodel\tstart\tstop\tprior_line\tassociation";
	static final String EVENT_HEADER="event_id\tinput_id\ttarget\tstart\tstop\tstrand\tmodel\tfamily\tmodel_from\tmodel_to\tmodel_length\tnative_truncation\tleft_edge\tright_edge\tprinted_bits\tGA_bits\tpercent_GA\tnative_included\tprinted_GA_disagreement\tmarginal_percent\tprinted_floor_pass\tfloor_uncertainty\ttier\treason\tsource_flags\tdecision_flags\tnative_table\tnative_line\tscore_basis\tpercent_precision";
	private static final int CONFIDENT=0, MARGINAL=1, INVALID=6;
	/** NcrnaReferenceTiers writes decision_flags with BRIDGE=16; never infer it from score. */
	private static final int BRIDGE_FLAG=16;
	private static final String[] TIERS={"CONFIDENT", "MARGINAL", "BELOW_MARGINAL", "UNCERTAIN_FLOOR", "REJECTED_WRONG_FAMILY", "COVERED_BY_PROTECTED_PRIOR", "INVALID_INTERVAL"};
	private static final int GROUP=0, PRIOR_OVERLAP=1, BRIDGE=2, WEAK=3;
	private static final String[] ASSOCIATIONS={"QUALIFYING_SAME_FAMILY_GROUP_MEMBER", "SAME_FAMILY_OVERLAP_NOT_RESCORING_PRIOR_BOUNDS", "TANDEM_BRIDGE_WITNESS_NOT_REPRESENTATIVE", "WEAK_OVERLAP_WITNESS_NOT_GROUPING_SUPPORT"};
}
