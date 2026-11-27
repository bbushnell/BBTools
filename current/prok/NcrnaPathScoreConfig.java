package prok;

import java.util.ArrayList;
import java.util.Arrays;

import map.ObjectIntMap;
import parse.LineParser1;
import prok.NcrnaPathGraphIO.Graph;
import prok.NcrnaTextIO.Input;

/** Explicit absolute coefficients for the saved-path cobalamin score experiment.
 * This reader does not fit, select or adopt coefficients. Each cell declares a
 * default for every registered member; optional exact-model rows override it.
 * Binding checks the actual loaded model axes before any graph is scored.
 * Float evaluation preserves the native A+B*L*identity*identity operation order;
 * the original c0/B1 score is retained bit-for-bit. The signed identity-minus-offset
 * term is squared without clamping; caller eligibility is unchanged. No labels or result objects
 * are allocated per candidate. After binding, coefficient arrays are immutable
 * and may be shared by caller workers. The explicitly selected runtime route
 * clamps the shifted identity at zero; saved historical graph evaluation above
 * retains its signed-square semantics.
 * @author Brian, Raiden
 */
public final class NcrnaPathScoreConfig{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaPathScoreConfig(final String path){
		final ObjectIntMap<String> byName=new ObjectIntMap<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, HEADER, path);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==5, "Score cells require cell, member, exact model or *, offset and absolute multiplier");
				final String name=NcrnaPathGraphIO.string(p, 0);
				int index=byName.get(name);
				if(index<0){
					index=cells.size();
					byName.put(name, index);
					cells.add(new Cell(name));
				}
				final int member=member(NcrnaPathGraphIO.string(p, 1));
				final String model=NcrnaPathGraphIO.string(p, 2);
				final float offset=p.parseFloat(3), multiplier=p.parseFloat(4);
				require((offset==0 || offset==.1f || offset==.2f || offset==.3f || offset==.4f || offset==.45f || offset==.5f || offset==.55f || offset==.6f)
					&& multiplier>0 && Float.isFinite(multiplier),
					"Declare one study offset 0/.1/.2/.3/.4/.45/.5/.55/.6 and a finite positive absolute multiplier; no relative reapplication");
				cells.get(index).rows.add(new Row(member, model, offset, multiplier));
			}
		}
		require(!cells.isEmpty(), "An evaluation must declare at least one explicit coefficient cell");
	}

	/** Bind once inside the verified native runtime callback, before its first graph.
	 * Defaults and overrides name actual model labels, never guessed model ordinals.
	 */
	public void bind(){
		bind(GeneCaller.ncrnaFamilies.toArray(new NcrnaFamily[GeneCaller.ncrnaFamilies.size()]));
	}
	/** Binds an explicit complete roster before a group is published to the caller. */
	public void bind(final NcrnaFamily[] families){
		require(!bound, "Each coefficient table binds exactly one native resource scope");
		require(families!=null, "Coefficient binding needs the actual loaded family roster");
		for(final NcrnaFamily family:families){
			require(family!=null, "A coefficient axis cannot contain an absent family");
			if(!NcrnaFamilyRegistry.isCobalamin(family.name)){continue;}
			final int m=member(family.name);
			require(names[m]==null && family.modelNames!=null, "Bind exactly one loaded model axis per target member");
			names[m]=family.modelNames;
			indices[m]=new ObjectIntMap<String>(String.class);
			for(int j=0; j<names[m].length; j++){
				require(!indices[m].contains(names[m][j]), "A target member cannot repeat an exact model label");
				indices[m].put(names[m][j], j);
			}
		}
		for(final String[] axis:names){require(axis!=null && axis.length>0, "All three target member axes must be loaded before scoring");}
		for(final Cell cell:cells){cell.bind();}
		bound=true;
	}
	public int cells(){return cells.size();}
	public Cell cell(final int index){return cells.get(index);}
	/** Select by the literal written cell name, never an assumed table ordinal. */
	public Cell cell(final String name){
		require(bound && name!=null && !name.isEmpty(), "Runtime selection needs a named cell on the actual bound model axes");
		for(final Cell cell:cells){if(cell.name.equals(name)){return cell;}}
		throw new IllegalArgumentException("Unknown coefficient cell: "+name);
	}

	/*--------------------------------------------------------------*/
	/*----------------             Cell             ----------------*/
	/*--------------------------------------------------------------*/

	public final class Cell{
		private Cell(final String name_){name=name_;}
		private void bind(){
			final boolean[] defaults=new boolean[MEMBERS.length];
			for(int m=0; m<MEMBERS.length; m++){
				offsets[m]=new float[names[m].length];
				multipliers[m]=new float[names[m].length];
			}
			for(final Row row:rows){
				if(!row.model.equals("*")){continue;}
				require(!defaults[row.member], "A cell must declare exactly one default per target member: "+name);
				defaults[row.member]=true;
				Arrays.fill(offsets[row.member], row.offset);
				Arrays.fill(multipliers[row.member], row.multiplier);
			}
			for(final boolean present:defaults){require(present, "Every cell must explicitly retain or change all three target members: "+name);}
			final boolean[][] overridden=new boolean[MEMBERS.length][];
			for(int m=0; m<MEMBERS.length; m++){overridden[m]=new boolean[names[m].length];}
			for(final Row row:rows){
				if(row.model.equals("*")){continue;}
				final int index=indices[row.member].get(row.model);
				require(index>=0 && !overridden[row.member][index], "Each model override must name one distinct actual model in its member: "+row.model);
				overridden[row.member][index]=true;
				offsets[row.member][index]=row.offset;
				multipliers[row.member][index]=row.multiplier;
			}
			rows.clear();
		}
		/** Assigns only target local scores after snapshot.restore(). Returns whether
		 * any score bits changed; a false result lets the driver reuse saved baseline
		 * statistics without a redundant native path search. The captured control
		 * must be A0/B1, the original coefficient scale for this study.
		 */
		public boolean apply(final Graph graph){
			require(bound && graph!=null && graph.snapshot!=null, "Score application requires the bound native axes and an affected restored graph");
			boolean changed=false;
			for(int i=0; i<graph.basisNodes.size; i++){
				final int m=member(graph.families.get(i));
				final int model=indices[m].get(graph.models.get(i));
				require(model>=0 && graph.additives.get(i)==0 && graph.multipliers.get(i)==1,
					"Offset coefficients are absolute on the bound original A0/B1 scale, not an already multiplied control");
				final float c=offsets[m][model], b=multipliers[m][model];
				if(c==0 && b==1){continue;}
				final float identity=graph.identities.get(i)-c;
				final float score=graph.additives.get(i)+b*graph.alignedLengths.get(i)*identity*identity;
				final long id=graph.basisNodes.get(i);
				if(Float.floatToRawIntBits(score)!=Float.floatToRawIntBits(graph.snapshot.originalScore(id))){
					graph.snapshot.setScore(id, score);
					changed=true;
				}
			}
			return changed;
		}
		/** Cold validation before attaching this immutable cell to a clan worker.
		 * Absolute study coefficients require the original A0/B1 member scale.
		 */
		void validateRuntime(final NcrnaFamily[] families){
			require(bound && families!=null && families.length==MEMBERS.length,
				"Runtime scoring requires all three bound cobalamin members");
			for(int m=0; m<MEMBERS.length; m++){
				final NcrnaFamily family=families[m];
				require(family!=null && MEMBERS[m].equals(family.name)
					&& Arrays.equals(names[m], family.modelNames) && family.scoreA==0 && family.scoreB==1,
					"Absolute runtime coefficients require the exact member/model order and original A0/B1 scale: "+MEMBERS[m]);
			}
		}
		/** Runtime-only formula after model selection, with no ranking or alignment.
		 * The owning clan validated its member/local model axes at installation.
		 */
		float runtimeScore(final int member, final int model, final int length, final float identity){
			assert(bound && member>=0 && member<MEMBERS.length && model>=0 && model<offsets[member].length)
				: shared.KillSwitch.assertDie("Clan scoring lost the member/local ordinal validated by validateRuntime: "+member+"/"+model);
			return clampedScore(offsets[member][model], multipliers[member][model], length, identity);
		}
		/** Immutable bound member-local offset, for an explicit saved-path search. */
		public float offset(final int member, final int model){
			checkModel(member, model);
			return offsets[member][model];
		}
		/** Immutable absolute coefficient, before a search's explicit relative factor. */
		public float multiplier(final int member, final int model){
			checkModel(member, model);
			return multipliers[member][model];
		}
		private void checkModel(final int member, final int model){
			if(!bound || member<0 || member>=MEMBERS.length || model<0 || model>=names[member].length){
				throw new IllegalArgumentException("Coefficient lookup needs the bound member/local model axis: "+member+"/"+model);
			}
		}
		public final String name;
		private final ArrayList<Row> rows=new ArrayList<Row>();
		private final float[][] offsets=new float[MEMBERS.length][], multipliers=new float[MEMBERS.length][];
	}
	/** Explicit runtime score on the original A0/B1 scale. Below-offset identity
	 * contributes zero instead of increasing again through a signed square.
	 */
	public static float clampedScore(final float offset, final float multiplier, final int length, final float identity){
		assert(length>0 && Float.isFinite(identity) && identity>=0 && identity<=1
			&& Float.isFinite(offset) && offset>=0 && offset<1 && multiplier>0 && Float.isFinite(multiplier))
			: shared.KillSwitch.assertDie("Runtime score operands lost the validated span/identity/absolute scale: length="+length
				+", identity="+identity+", offset="+offset+", multiplier="+multiplier);
		final float shifted=Math.max(0f, identity-offset);
		final float score=0f+multiplier*length*shifted*shifted;
		assert(Float.isFinite(score)) : shared.KillSwitch.assertDie("Runtime path scoring must remain finite on caller workers: length="+length
			+", identity="+identity+", offset="+offset+", multiplier="+multiplier);
		return score;
	}
	private static final class Row{
		Row(final int member_, final String model_, final float offset_, final float multiplier_){
			member=member_;
			model=model_;
			offset=offset_;
			multiplier=multiplier_;
		}
		final int member;
		final String model;
		final float offset, multiplier;
	}
	private static int member(final String family){
		for(int i=0; i<MEMBERS.length; i++){if(MEMBERS[i].equals(family)){return i;}}
		throw new IllegalArgumentException("This score experiment requires an explicitly registered cobalamin member: "+family);
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final ArrayList<Cell> cells=new ArrayList<Cell>();
	private final String[][] names=new String[MEMBERS.length][];
	@SuppressWarnings("unchecked")
	private final ObjectIntMap<String>[] indices=new ObjectIntMap[MEMBERS.length];
	private boolean bound;
	public static final String HEADER="cell\tmember\tmodel\tidentity_offset\tabsolute_multiplier";
	private static final String[] MEMBERS={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
}
