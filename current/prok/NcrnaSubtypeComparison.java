package prok;

import java.util.ArrayList;
import java.util.Arrays;

import map.IntHashMap;
import map.IntLongHashMap;
import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.IntList;

/** Optional saved cobalamin subtype comparisons for SummarizeNcrnaBenchmark.
 * Uses exact native input/parent keys and original query-line witnesses, never
 * overlap matching. Subtype and union assignments are independent: a subtype
 * loss may retain a union match to another subtype OR to the same subtype.
 * No caller, regrading, model choice or biological subtype inference occurs.
 * Single-threaded: parser rows are borrowed only during each method call;
 * retained identity keys, model names and geometry are owned copies. Mutable
 * lookup probes and output buffers are reused and never stored as map keys.
 * @author Ganyu
 */
final class NcrnaSubtypeComparison{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	NcrnaSubtypeComparison(final int count){
		require(count>0, "The benchmark must declare at least one observation cell");
		cells=new Cell[count];
		for(int i=0; i<count; i++){cells[i]=new Cell();}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Native readers        ----------------*/
	/*--------------------------------------------------------------*/

	/** Feed the already parsed ordinary raw report, including excluded queries. */
	void ordinary(final int cell, final LineParser1 p, final int offset){
		if(!p.termEquals("cobalamin", offset+2)){return;}
		final Cell c=cells[cell];
		if(p.termEquals("REFERENCE", offset+1)){
			final int id=parent(p, offset, offset+8, true);
			final Parent parent=parents.get(id);
			final Target shape=geometry(p, offset+3);
			if(parent.geometry==null){parent.geometry=shape.copy();}
			else{require(parent.geometry.equals(shape), "Union reference geometry changed for an exact input/parent key");}
			require(!c.states.contains(id), "Duplicate union reference parent in one cell");
			final int state=referenceState(p, offset+7);
			c.states.put(id, state);
			c.parentIds.add(id);
			if(state==TP){c.matches.put(id, positive(p, offset+11));}
			else{require(p.termEquals(".", offset+11), "Unmatched union reference must not name a query line");}
		}else{
			require(p.termEquals("QUERY", offset+1), "Ordinary rows must be QUERY or REFERENCE");
			final long line=positive(p, offset+9);
			final Target key=queryKey(p, offset, line);
			require(c.queryIndex.get(key)<0, "Original query line repeats within one input/cell");
			final Target retained=key.copy();
			final int model=model(p, offset+15);
			final boolean included=!p.termEquals("EXCLUDED", offset+7);
			final int matched=p.termEquals("TP", offset+7) ? parent(p, offset, offset+10, true) : -1;
			final Query query=new Query(model, geometry(p, offset+3).copy(), line, matched, included);
			c.queryIndex.put(retained, c.queries.size());
			c.queries.add(query);
		}
	}

	/** Native subtype totals prevent a truncated detail table from looking stable. */
	void summary(final int cell, final LineParser1 p, final int offset, final int inputs, final double overlap){
		final View v=cells[cell].views[slot(p, offset)];
		require(v.expected==null, "Each cell needs one summary per cobalamin subtype and view");
		require(p.parseLong(offset+2)==inputs, "Subtype and ordinary input populations disagree");
		v.expected=new long[COUNTS];
		for(int i=0; i<COUNTS; i++){
			v.expected[i]=p.parseLong(offset+3+i);
			require(v.expected[i]>=0, "Subtype observation counts must be nonnegative");
		}
		require(Math.addExact(v.expected[TP_COUNT], v.expected[FN_COUNT])==v.expected[SCORED],
			"Subtype TP plus FN must equal its scored parent count");
		require(v.expected[BOTH]<=v.expected[EXACT5] && v.expected[BOTH]<=v.expected[EXACT3]
			&& v.expected[EXACT5]<=v.expected[TP_COUNT] && v.expected[EXACT3]<=v.expected[TP_COUNT],
			"Subtype exact endpoints must refer to matched parents");
		ratio(p, offset+13, v.expected[TP_COUNT], Math.addExact(v.expected[TP_COUNT], v.expected[FP_COUNT]));
		ratio(p, offset+14, v.expected[TP_COUNT], v.expected[SCORED]);
		require(Double.compare(p.parseDouble(offset+15), overlap)==0,
			"Subtype and union reports must use the same overlap rule");
	}

	/** Read both roles so every subtype TP retains its actual original query. */
	void raw(final int cell, final LineParser1 p, final int offset){
		final int slot=slot(p, offset);
		final Cell c=cells[cell];
		final View v=c.views[slot];
		if(p.termEquals("REFERENCE", offset+3)){
			final int id=parent(p, offset+2, offset+4, false);
			require(c.states.contains(id), "Subtype parent is absent from this cell's ordinary union");
			require(parents.get(id).geometry.equals(geometry(p, offset+7)),
				"Subtype and union reference coordinates disagree for the same parent");
			require(!v.states.contains(id), "Subtype reference parent repeats in one cell/view");
			final int state=referenceState(p, offset+11);
			v.states.put(id, state);
			v.parentIds.add(id);
			if(state==MASK){v.observed[MASKED_REF]++;}
			else{
				v.observed[SCORED]++;
				v.observed[state==TP ? TP_COUNT : FN_COUNT]++;
			}
			if(state==TP){v.references.put(id, positive(p, offset+13));}
			else{require(p.termEquals(".", offset+13), "Unmatched subtype parent must not name a query line");}
		}else{
			require(p.termEquals("QUERY", offset+3), "Subtype rows must be QUERY or REFERENCE");
			final long line=positive(p, offset+5);
			final int qi=c.queryIndex.get(queryKey(p, offset+2, line));
			require(qi>=0, "Subtype query line has no original ordinary query witness");
			final Query q=c.queries.get(qi);
			require(q.included && models.get(q.model).type==slot%TYPES.length,
				"Subtype query must be an included original query of this model subtype");
			probe.set(p, offset+6);
			require(models.get(q.model).name.equals(probe) && q.geometry.equals(geometry(p, offset+7)),
				"Subtype query model or coordinates differ from the original union query");
			require(!v.seenQueries.contains(qi), "Original query repeats within one subtype view");
			v.seenQueries.put(qi, 1);
			if(p.termEquals("true", offset+15)){v.observed[DUPLICATES]++;}
			else{require(p.termEquals("false", offset+15), "Subtype query duplicate flag must be explicit");}
			if(p.termEquals("TP", offset+11)){
				final int id=parent(p, offset+2, offset+12, false);
				require(!v.queries.contains(id), "Subtype queries consume a parent more than once");
				v.queries.put(id, line);
				final boolean exact5=p.parseInt(offset+20)==0;
				final boolean exact3=p.parseInt(offset+21)==0;
				if(exact5){v.observed[EXACT5]++;}
				if(exact3){v.observed[EXACT3]++;}
				if(exact5 && exact3){v.observed[BOTH]++;}
			}else{
				require(p.termEquals(".", offset+12), "Non-TP subtype query must not name a scored match");
				if(p.termEquals("FP", offset+11)){v.observed[FP_COUNT]++;}
				else{
					require(p.termEquals("MASKED", offset+11), "Unknown subtype query outcome");
					v.observed[MASKED_CALL]++;
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Validation           ----------------*/
	/*--------------------------------------------------------------*/

	void validate(){
		for(int cell=0; cell<cells.length; cell++){
			final Cell c=cells[cell];
			for(int i=0; i<c.parentIds.size; i++){
				final int id=c.parentIds.get(i);
				if(c.states.get(id)==TP){
					final Query q=unionQuery(c, id);
					require(q.parent==id && q.included, "Union reference and query must name the same exact parent");
				}
			}
			for(final Query q:c.queries){
				if(q.parent>=0){
					require(c.states.get(q.parent)==TP && c.matches.get(q.parent)==q.line,
						"Union TP query must have a reciprocal reference witness");
				}
			}
			for(int slot=0; slot<c.views.length; slot++){
				final View v=c.views[slot];
				require(v.expected!=null, "A declared cell is missing a cobalamin subtype/view summary");
				if(!Arrays.equals(v.expected, v.observed)){
					throw new IllegalArgumentException("Subtype raw rows disagree with native summary counts: cell="+cell+
						" view="+NcrnaSubtypeGrade.VIEWS[slot/TYPES.length]+" subtype="+TYPES[slot%TYPES.length]+
						" expected="+Arrays.toString(v.expected)+" observed="+Arrays.toString(v.observed));
				}
				require(v.references.size()==v.queries.size(), "Subtype TP reference/query populations differ");
				v.parentIds.sort();
				for(int i=0; i<v.parentIds.size; i++){
					final int id=v.parentIds.get(i);
					if(v.states.get(id)==TP){
						require(v.queries.contains(id) && v.queries.get(id)==v.references.get(id),
							"Subtype TP reference/query line witnesses disagree");
					}
				}
				for(int i=0; i<c.queries.size(); i++){
					final Query q=c.queries.get(i);
					if(q.included && models.get(q.model).type==slot%TYPES.length){
						require(v.seenQueries.contains(i), "Subtype view omitted an original query of its model type");
					}
				}
			}
		}
		for(final Cell c:cells){
			for(int slot=0; slot<c.views.length; slot++){
				final View a=cells[0].views[slot];
				final View b=c.views[slot];
				require(a.parentIds.size==b.parentIds.size, "Compared subtype views have different parent universes");
				for(int i=0; i<a.parentIds.size; i++){
					final int id=a.parentIds.get(i);
					require(id==b.parentIds.get(i) && (a.states.get(id)==MASK)==(b.states.get(id)==MASK),
						"Compared subtype views must preserve exact parents and scored/masked status");
				}
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Exact comparisons       ----------------*/
	/*--------------------------------------------------------------*/

	/** LOST assesses the right union; GAINED assesses the left union. Both witnesses
	 * are emitted so a view transition is never silently called a biological loss.
	 */
	void compare(final ByteBuilder context, final int left, final int right, final Output summary, final Output changes){
		assert(summary!=null && changes!=null) : "Optional comparison writers share the ordinary report's checked output lifetime";
		final Cell a=cells[left];
		final Cell b=cells[right];
		for(int slot=0; slot<a.views.length; slot++){
			Arrays.fill(counts, 0);
			final View av=a.views[slot];
			final View bv=b.views[slot];
			for(int i=0; i<av.parentIds.size; i++){
				final int id=av.parentIds.get(i);
				final int before=av.states.get(id);
				final int after=bv.states.get(id);
				if((before==TP)==(after==TP)){continue;}
				final boolean lost=before==TP;
				final int category=category(lost ? b : a, id, slot%TYPES.length);
				counts[(lost ? 0 : CATEGORIES.length)+category]++;
				start(context, slot).tab().append(lost ? "LOST" : "GAINED").tab()
					.append(parents.get(id).key.bytes).tab().append(parents.get(id).geometry.bytes).tab()
					.append(STATES[before]).tab().append(STATES[after]).tab();
				witness(a, id).tab();
				witness(b, id).tab().append(lost ? "right" : "left").tab().append(CATEGORIES[category]).nl();
				changes.writer.print(buffer);
			}
			long lost=0;
			long gained=0;
			for(int i=0; i<CATEGORIES.length; i++){
				lost+=counts[i];
				gained+=counts[CATEGORIES.length+i];
			}
			require(av.expected[TP_COUNT]-lost+gained==bv.expected[TP_COUNT], "Exact subtype set changes must reconcile with both native TP totals");
			start(context, slot).tab().append(av.expected[TP_COUNT]).tab().append(bv.expected[TP_COUNT])
				.tab().append(lost).tab().append(gained);
			for(final long count:counts){buffer.tab().append(count);}
			summary.writer.print(buffer.append("\tSAVED_MATCH_WITNESSES_NONADDITIVE_NOT_BIOLOGICAL_SUBTYPE_PROOF\n"));
		}
	}

	private ByteBuilder start(final ByteBuilder context, final int slot){
		assert(slot>=0 && slot<SLOTS) : "Validated view/subtype pairs index fixed native cobalamin views";
		return buffer.clear().append(context).tab().append(NcrnaSubtypeGrade.VIEWS[slot/TYPES.length]).tab().append(TYPES[slot%TYPES.length]);
	}
	private ByteBuilder witness(final Cell c, final int id){
		final int state=c.states.get(id);
		require(state>0, "Every subtype parent needs an ordinary union outcome");
		buffer.append(STATES[state]).tab();
		if(state==TP){
			final Query q=unionQuery(c, id);
			buffer.append(q.line).tab().append(models.get(q.model).name.bytes);
		}else{buffer.append(".\t.");}
		return buffer;
	}
	private int category(final Cell c, final int id, final int subtype){
		final int state=c.states.get(id);
		if(state==FN){return 2;}
		if(state==MASK){return 3;}
		require(state==TP, "An ordinary union outcome is required to classify a view change");
		return models.get(unionQuery(c, id).model).type==subtype ? 1 : 0;
	}
	private Query unionQuery(final Cell c, final int id){
		final Parent parent=parents.get(id);
		final long line=c.matches.get(id);
		require(line>0, "A union TP requires an original matched query line");
		key.clear().append(parent.key.bytes, 0, parent.inputLength).tab().append(line);
		probe.set(key.array, 0, key.length());
		final int index=c.queryIndex.get(probe);
		require(index>=0, "Union matched query line is missing from the saved raw report");
		return c.queries.get(index);
	}
	/*--------------------------------------------------------------*/
	/*----------------       Keys and parsing       ----------------*/
	/*--------------------------------------------------------------*/

	private int parent(final LineParser1 p, final int input, final int field, final boolean create){
		require(!p.termEquals(".", field), "A compared parent needs its stable native identifier");
		key.clear();
		append(key, p, input);
		final int inputLength=key.length();
		require(inputLength>0, "Parent input identity must be nonempty");
		key.tab();
		append(key, p, field);
		require(key.length()>inputLength+1, "Parent identity must be nonempty");
		probe.set(key.array, 0, key.length());
		int id=parentIndex.get(probe);
		if(id<0){
			require(create, "Subtype parent is absent from the ordinary raw report");
			final Target retained=probe.copy();
			id=parents.size();
			parents.add(new Parent(retained, inputLength));
			parentIndex.put(retained, id);
		}
		return id;
	}
	private Target queryKey(final LineParser1 p, final int input, final long line){
		assert(line>0) : "Only positive physical GFF query-line witnesses form lookup keys";
		key.clear();
		append(key, p, input).tab().append(line);
		probe.set(key.array, 0, key.length());
		return probe;
	}
	private Target geometry(final LineParser1 p, final int seqid){
		key.clear();
		for(int i=0; i<4; i++){
			if(i>0){key.tab();}
			append(key, p, seqid+i);
		}
		probe.set(key.array, 0, key.length());
		return probe;
	}
	private int model(final LineParser1 p, final int field){
		probe.set(p, field);
		int id=modelIndex.get(probe);
		if(id>=0){return id;}
		int type=-1;
		for(int i=0; i<PREFIXES.length; i++){if(p.termStartsWith(PREFIXES[i], field)){type=i;}}
		require(type>=0, "Unknown ordinary cobalamin model prefix cannot supply subtype evidence");
		final Target retained=probe.copy();
		id=models.size();
		models.add(new Model(retained, type));
		modelIndex.put(retained, id);
		return id;
	}
	private static ByteBuilder append(final ByteBuilder out, final LineParser1 p, final int field){
		p.setBounds(field);
		return out.append(p.line(), p.a(), p.b()-p.a());
	}
	private static int slot(final LineParser1 p, final int offset){
		int view=-1;
		int type=-1;
		for(int i=0; i<NcrnaSubtypeGrade.VIEWS.length; i++){if(p.termEquals(NcrnaSubtypeGrade.VIEWS[i], offset)){view=i;}}
		for(int i=0; i<TYPES.length; i++){if(p.termEquals(TYPES[i], offset+1)){type=i;}}
		require(view>=0 && type>=0, "Unknown native cobalamin subtype or view");
		return view*TYPES.length+type;
	}
	private static int referenceState(final LineParser1 p, final int field){
		if(p.termEquals("TP", field)){return TP;}
		if(p.termEquals("FN", field)){return FN;}
		require(p.termEquals("DONT_CARE", field), "Unknown saved reference outcome");
		return MASK;
	}
	private static long positive(final LineParser1 p, final int field){
		final long value=p.parseLong(field);
		require(value>0, "Matched query lines must be positive physical GFF lines");
		return value;
	}
	private static void ratio(final LineParser1 p, final int field, final long numerator, final long denominator){
		if(denominator==0){require(p.termEquals("NA", field), "A zero subtype ratio denominator must be NA");}
		else{
			final double value=p.parseDouble(field);
			require(Double.isFinite(value) && Math.abs(value-numerator/(double)denominator)<=0.000000500001,
				"Subtype ratio disagrees with native counts at six-place precision");
		}
	}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------      Owned observations      ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Cell{
		Cell(){for(int i=0; i<views.length; i++){views[i]=new View();}}
		final IntHashMap states=new IntHashMap();
		final IntLongHashMap matches=new IntLongHashMap();
		final IntList parentIds=new IntList();
		final ObjectIntMap<Target> queryIndex=new ObjectIntMap<Target>(Target.class);
		final ArrayList<Query> queries=new ArrayList<Query>();
		final View[] views=new View[SLOTS];
	}
	private static final class View{
		final IntHashMap states=new IntHashMap();
		final IntHashMap seenQueries=new IntHashMap();
		final IntLongHashMap references=new IntLongHashMap();
		final IntLongHashMap queries=new IntLongHashMap();
		final IntList parentIds=new IntList();
		final long[] observed=new long[COUNTS];
		long[] expected;
	}
	private static final class Parent{
		Parent(final Target key_, final int inputLength_){
			key=key_;
			inputLength=inputLength_;
		}
		final Target key;
		final int inputLength;
		Target geometry;
	}
	private static final class Query{
		Query(final int model_, final Target geometry_, final long line_, final int parent_, final boolean included_){
			model=model_;
			geometry=geometry_;
			line=line_;
			parent=parent_;
			included=included_;
		}
		final Target geometry;
		final int model, parent;
		final long line;
		final boolean included;
	}
	private static final class Model{
		Model(final Target name_, final int type_){
			name=name_;
			type=type_;
		}
		final Target name;
		final int type;
	}
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Cell[] cells;
	private final ObjectIntMap<Target> parentIndex=new ObjectIntMap<Target>(Target.class);
	private final ObjectIntMap<Target> modelIndex=new ObjectIntMap<Target>(Target.class);
	private final ArrayList<Parent> parents=new ArrayList<Parent>();
	private final ArrayList<Model> models=new ArrayList<Model>();
	private final ByteBuilder key=new ByteBuilder();
	private final ByteBuilder buffer=new ByteBuilder();
	private final Target probe=new Target();
	private final long[] counts=new long[8];

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final int TP=1, FN=2, MASK=3;
	private static final String[] STATES={"MISSING", "TP", "FN", "DONT_CARE"};
	private static final String[] TYPES={"RF00174", "RF01482", "RF01689"};
	private static final String[] PREFIXES={"cobalamin_RF00174", "cobalamin_RF01482", "cobalamin_RF01689"};
	private static final String[] CATEGORIES={"OTHER_SUBTYPE_MODEL", "SAME_SUBTYPE_REMATCH", "UNION_MISS", "UNION_MASKED"};
	private static final int SLOTS=6, COUNTS=10;
	private static final int TP_COUNT=0, FP_COUNT=1, FN_COUNT=2, SCORED=3, MASKED_REF=4, MASKED_CALL=5,
		DUPLICATES=6, EXACT5=7, EXACT3=8, BOTH=9;
	static final String SUMMARY_TAIL="view\tsubtype\tleft_TP\tright_TP\tlost_parents\tgained_parents\tlost_other_model\tlost_same_model\tlost_union_miss\tlost_union_masked\tgained_other_model\tgained_same_model\tgained_union_miss\tgained_union_masked\tqualification";
	static final String CHANGE_TAIL="view\tsubtype\tchange\tinput_id\tparent_id\tseqid\tstart\tstop\tstrand\tleft_subtype_outcome\tright_subtype_outcome\tleft_union_outcome\tleft_union_query_line\tleft_union_model\tright_union_outcome\tright_union_query_line\tright_union_model\tassessed_side\tunion_class";
}
