package test.tax;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicLong;

import structures.IntList;
import tax.GiToTaxid;
import tax.TaxNode;
import tax.TaxServer;
import tax.TaxTree;

/** Finite oracle for TaxServer missing-query accounting.
 * @author Brian Bushnell, Jean */
public final class TaxServerNotFoundAccountingTest{

	/** Prevents instantiation. */
	private TaxServerNotFoundAccountingTest(){}

	/** Runs baseline or fixed expectations for missing-query accounting.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If reflective setup or method calls fail */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new AssertionError("Unknown argument: "+arg);}
		}

		final HashMap<String, Integer> typeMap=makeTypeMap();
		final int gi=requireAlias(typeMap, "gi");
		final int name=requireAlias(typeMap, "name");
		final int taxid=requireAlias(typeMap, "taxid");
		final TaxTree tree=tinyTree();
		final AtomicLong notFound=new AtomicLong(0);
		final TaxServer server=server(tree, notFound);
		final Method toIntList=TaxServer.class.getDeclaredMethod("toIntList",
				int.class, String[].class);
		toIntList.setAccessible(true);
		final Method getTaxNodeGi=TaxServer.class.getDeclaredMethod("getTaxNodeGi",
				long.class);
		getTaxNodeGi.setAccessible(true);
		final Method getTaxNodeByName=TaxServer.class.getDeclaredMethod(
				"getTaxNodeByName", String.class);
		getTaxNodeByName.setAccessible(true);
		final Method getTaxNodeTaxid=TaxServer.class.getDeclaredMethod(
				"getTaxNodeTaxid", int.class);
		getTaxNodeTaxid.setAccessible(true);

		try(GiState ignored=new GiState()){
			requireDirectHit(getTaxNodeGi, server, notFound, Long.valueOf(HIT_GI));
			requireDirectHit(getTaxNodeByName, server, notFound, HIT_NAME);
			requireDirectHit(getTaxNodeTaxid, server, notFound, Integer.valueOf(HIT_TAXID));

			requireDirectMiss(getTaxNodeGi, server, notFound, Long.valueOf(MISS_GI));
			requireDirectMiss(getTaxNodeGi, server, notFound,
					Long.valueOf(MAPPED_MISS_GI));
			requireDirectMiss(getTaxNodeByName, server, notFound, MISS_NAME);
			requireDirectMiss(getTaxNodeTaxid, server, notFound, Integer.valueOf(MISS_TAXID));

			requireListHit(toIntList, server, notFound, gi, Long.toString(HIT_GI));
			requireListHit(toIntList, server, notFound, name, HIT_NAME);
			requireListHit(toIntList, server, notFound, taxid, Integer.toString(HIT_TAXID));

			final long expectedMissing=baseline ? 2 : 1;
			requireListMiss(toIntList, server, notFound, gi, Long.toString(MISS_GI), expectedMissing);
			requireListMiss(toIntList, server, notFound, gi,
					Long.toString(MAPPED_MISS_GI), expectedMissing);
			requireListMiss(toIntList, server, notFound, name, MISS_NAME, expectedMissing);
			requireListMiss(toIntList, server, notFound, taxid,
					Integer.toString(MISS_TAXID), expectedMissing);
		}

		System.out.println("TAXSERVER_NOTFOUND_ACCOUNTING_TEST_PASS expect="+
				(baseline ? "baseline" : "fixed"));
	}

	/** Reflects TaxServer.makeTypeMap.
	 * @return Route alias map
	 * @throws Exception If the method cannot be reflected or invoked */
	@SuppressWarnings("unchecked")
	private static HashMap<String, Integer> makeTypeMap() throws Exception{
		final Method method=TaxServer.class.getDeclaredMethod("makeTypeMap");
		method.setAccessible(true);
		return (HashMap<String, Integer>)method.invoke(null);
	}

	/** Allocates a TaxServer without running the file-loading constructor.
	 * @param tree Tiny explicit tree
	 * @param notFound Missing-query counter
	 * @return Initialized server with only fields needed by this fixture
	 * @throws Exception If Unsafe allocation or reflection fails */
	private static TaxServer server(TaxTree tree, AtomicLong notFound) throws Exception{
		Class<?> unsafeClass=Class.forName("sun.misc.Unsafe");
		Field unsafeField=unsafeClass.getDeclaredField("theUnsafe");
		unsafeField.setAccessible(true);
		Object unsafe=unsafeField.get(null);
		Method allocateInstance=unsafeClass.getMethod("allocateInstance", Class.class);

		TaxServer server=(TaxServer)allocateInstance.invoke(unsafe, TaxServer.class);
		setField(TaxServer.class, server, "tree", tree);
		setField(TaxServer.class, server, "commonMap", new HashMap<String, String>());
		setField(TaxServer.class, server, "notFound", notFound);
		return server;
	}

	/** Sets one reflected field.
	 * @param owner Class declaring the field
	 * @param target Instance to modify, or null for static fields
	 * @param name Field name
	 * @param value New field value
	 * @throws Exception If the field cannot be reflected or set */
	private static void setField(Class<?> owner, Object target, String name,
			Object value) throws Exception{
		final Field field=owner.getDeclaredField(name);
		field.setAccessible(true);
		field.set(target, value);
	}

	/** Reads one reflected field.
	 * @param owner Class declaring the field
	 * @param target Instance to read, or null for static fields
	 * @param name Field name
	 * @return Current field value
	 * @throws Exception If the field cannot be reflected or read */
	private static Object getField(Class<?> owner, Object target, String name)
			throws Exception{
		final Field field=owner.getDeclaredField(name);
		field.setAccessible(true);
		return field.get(target);
	}

	/** Builds the smallest tree needed by TaxServer GI, NAME and TAXID lookup.
	 * @return Tiny in-memory taxonomy tree with one resolvable child node */
	private static TaxTree tinyTree(){
		TaxNode[] nodes=new TaxNode[MAPPED_MISS_TAXID+1];
		nodes[1]=new TaxNode(1, 1, TaxTree.LIFE, TaxTree.LIFE_E, "Life");
		nodes[HIT_TAXID]=new TaxNode(HIT_TAXID, 1, TaxTree.SPECIES,
				TaxTree.SPECIES_E, HIT_NAME);
		TaxTree tree=new TaxTree(nodes, null, 0, false, false, false,
				TaxTree.NO_RANK_E);
		tree.hashNames(false);
		return tree;
	}

	/** Requires one direct helper call to resolve without counter drift.
	 * @param method Reflected TaxServer helper
	 * @param server Initialized server
	 * @param notFound Missing-query counter
	 * @param query Query argument
	 * @throws Exception If the helper cannot be invoked */
	private static void requireDirectHit(Method method, TaxServer server,
			AtomicLong notFound, Object query) throws Exception{
		final long before=notFound.get();
		final TaxNode tn=(TaxNode)method.invoke(server, query);
		check(tn!=null, "Expected direct "+method.getName()+" hit for "+query);
		check(tn.id==HIT_TAXID,
				method.getName()+" returned ID "+tn.id+" for "+query);
		check(notFound.get()==before,
				method.getName()+" hit changed notFound from "+before+
				" to "+notFound.get());
	}

	/** Requires one direct helper call to miss and increment exactly once.
	 * @param method Reflected TaxServer helper
	 * @param server Initialized server
	 * @param notFound Missing-query counter
	 * @param query Query argument
	 * @throws Exception If the helper cannot be invoked */
	private static void requireDirectMiss(Method method, TaxServer server,
			AtomicLong notFound, Object query) throws Exception{
		final long before=notFound.get();
		final TaxNode tn=(TaxNode)method.invoke(server, query);
		check(tn==null, "Expected direct "+method.getName()+" miss for "+query);
		check(notFound.get()==before+1,
				method.getName()+" miss should increment once; before="+
				before+", after="+notFound.get());
	}

	/** Requires TaxServer.toIntList to resolve one ID without counter drift.
	 * @param method Reflected toIntList method
	 * @param server Initialized server
	 * @param notFound Missing-query counter
	 * @param type Encoded route type
	 * @param name Query term
	 * @throws Exception If toIntList cannot be invoked */
	private static void requireListHit(Method method, TaxServer server,
			AtomicLong notFound, int type, String name) throws Exception{
		final long before=notFound.get();
		requireList(toIntList(method, server, type, name), HIT_TAXID);
		check(notFound.get()==before,
				"toIntList hit changed notFound from "+before+
				" to "+notFound.get());
	}

	/** Requires TaxServer.toIntList to miss and increment the expected amount.
	 * @param method Reflected toIntList method
	 * @param server Initialized server
	 * @param notFound Missing-query counter
	 * @param type Encoded route type
	 * @param name Query term
	 * @param expected Expected notFound delta
	 * @throws Exception If toIntList cannot be invoked */
	private static void requireListMiss(Method method, TaxServer server,
			AtomicLong notFound, int type, String name, long expected)
			throws Exception{
		final long before=notFound.get();
		requireList(toIntList(method, server, type, name));
		final long delta=notFound.get()-before;
		check(delta==expected,
				"toIntList miss should increment "+expected+" time(s); found "+delta);
	}

	/** Runs TaxServer.toIntList reflectively.
	 * @param method Reflected toIntList method
	 * @param server Initialized server
	 * @param type Encoded route type
	 * @param names Query names
	 * @return Resolved TaxIDs
	 * @throws Exception If toIntList cannot be invoked or throws */
	private static IntList toIntList(Method method, TaxServer server, int type,
			String... names) throws Exception{
		try{
			return (IntList)method.invoke(server, type, (Object)names);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			if(cause instanceof Error){throw (Error)cause;}
			if(cause instanceof Exception){throw (Exception)cause;}
			throw new RuntimeException(cause);
		}
	}

	/** Requires an alias to exist.
	 * @param typeMap TaxServer alias map
	 * @param alias Route alias to inspect
	 * @return Actual route value */
	private static int requireAlias(final HashMap<String, Integer> typeMap,
			final String alias){
		final Integer actual=typeMap.get(alias);
		check(actual!=null, "Missing route alias "+alias);
		return actual.intValue();
	}

	/** Requires an IntList to contain exactly the expected IDs.
	 * @param list List to inspect
	 * @param expected Expected IDs */
	private static void requireList(IntList list, int... expected){
		check(list.size==expected.length,
				"Expected "+expected.length+" IDs; found "+list.size);
		for(int i=0; i<expected.length; i++){
			check(list.get(i)==expected[i],
					"Expected ID "+expected[i]+" at "+i+"; found "+list.get(i));
		}
	}

	/** Test checks stay active under -da.
	 * @param condition Condition required for test success
	 * @param message Diagnostic message for failure */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Saves, replaces and restores global GI lookup state. */
	private static final class GiState implements AutoCloseable{

		/** Installs a tiny explicit GI table.
		 * @throws Exception If GI static state cannot be reflected */
		GiState() throws Exception{
			oldDiskMode=GiToTaxid.DISK_MODE;
			oldArray=(int[][])getField(GiToTaxid.class, null, "array");
			oldMaxGiLoaded=((Long)getField(GiToTaxid.class, null,
					"maxGiLoaded")).longValue();
			oldDiskTable=getField(GiToTaxid.class, null, "diskTable");
			oldFileString=(String)getField(GiToTaxid.class, null, "fileString");
			oldInitialized=((Boolean)getField(GiToTaxid.class, null,
					"initialized")).booleanValue();

			int[][] array=new int[1][];
			array[0]=new int[MAPPED_MISS_GI+1];
			array[0][HIT_GI]=HIT_TAXID;
			array[0][MAPPED_MISS_GI]=MAPPED_MISS_TAXID;
			setField(GiToTaxid.class, null, "diskTable", null);
			setField(GiToTaxid.class, null, "array", array);
			setField(GiToTaxid.class, null, "maxGiLoaded",
					Long.valueOf(MAPPED_MISS_GI));
			setField(GiToTaxid.class, null, "fileString", getClass().getName());
			setField(GiToTaxid.class, null, "DISK_MODE", Boolean.FALSE);
			setField(GiToTaxid.class, null, "initialized", Boolean.TRUE);
		}

		/** Restores the exact GI state replaced by this fixture.
		 * @throws RuntimeException If reflection fails */
		@Override
		public void close(){
			try{
				setField(GiToTaxid.class, null, "diskTable", oldDiskTable);
				setField(GiToTaxid.class, null, "array", oldArray);
				setField(GiToTaxid.class, null, "maxGiLoaded",
						Long.valueOf(oldMaxGiLoaded));
				setField(GiToTaxid.class, null, "fileString", oldFileString);
				setField(GiToTaxid.class, null, "DISK_MODE", Boolean.valueOf(oldDiskMode));
				setField(GiToTaxid.class, null, "initialized",
						Boolean.valueOf(oldInitialized));
			}catch(Exception e){
				throw new RuntimeException("Failed to restore GI state.", e);
			}
		}

		/** Original disk-mode flag. */
		private final boolean oldDiskMode;
		/** Original in-memory GI table. */
		private final int[][] oldArray;
		/** Original highest loaded GI. */
		private final long oldMaxGiLoaded;
		/** Original disk table object. */
		private final Object oldDiskTable;
		/** Original table source descriptor. */
		private final String oldFileString;
		/** Original initialized flag. */
		private final boolean oldInitialized;
	}

	/** GI that resolves through the tiny explicit table. */
	private static final int HIT_GI=7;
	/** GI that returns 0 from the tiny explicit table. */
	private static final int MISS_GI=8;
	/** GI that maps to an absent positive TaxID. */
	private static final int MAPPED_MISS_GI=9;
	/** TaxID to resolve. */
	private static final int HIT_TAXID=2;
	/** TaxID absent from the tiny tree. */
	private static final int MISS_TAXID=0;
	/** Positive TaxID absent from the tiny tree. */
	private static final int MAPPED_MISS_TAXID=3;
	/** Taxon name to resolve. */
	private static final String HIT_NAME="Tiny species";
	/** Taxon name absent from the tiny tree and common map. */
	private static final String MISS_NAME="missing species";
}
