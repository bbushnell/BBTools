package test.tax;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicLong;

import structures.IntList;
import tax.TaxNode;
import tax.TaxServer;
import tax.TaxTree;

/** Finite oracle for TaxServer header ancestor conversion.
 * @author Brian Bushnell, Jean */
public final class TaxServerHeaderAncestorTest{

	/** Runs baseline or fixed expectations for header ancestor conversion.
	 * @param args Optional expect=baseline or expect=fixed
	 * @throws Exception If reflective setup or method calls fail */
	public static void main(String[] args) throws Exception{
		boolean baseline=false;
		for(String arg : args){
			if(arg.equals("expect=baseline")){baseline=true;}
			else if(arg.equals("expect=fixed")){baseline=false;}
			else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}

		final HashMap<String, Integer> typeMap=typeMap();
		final int taxid=requireAlias(typeMap, "id", EXPECTED_TAXID);
		final int header=requireAlias(typeMap, "header", EXPECTED_HEADER);
		final int silva=requireAlias(typeMap, "silvaheader", EXPECTED_SILVAHEADER);
		final int ptHeader=requireAlias(typeMap, "pt_header", EXPECTED_PT_HEADER);
		final int scSilva=requireAlias(typeMap, "sc_silvaheader", EXPECTED_SC_SILVAHEADER);

		final Method toIntList=TaxServer.class.getDeclaredMethod("toIntList",
				int.class, String[].class);
		toIntList.setAccessible(true);

		final AtomicLong notFound=new AtomicLong(0);
		final TaxServer server=server(tinyTree(), notFound);
		requireList(toIntList(toIntList, server, taxid, "2", "2"), 2, 2);
		observeMissingTaxid(toIntList, server, taxid, notFound);

		if(baseline){
			requireBadType(toIntList, server, header, "tid|2|ncbi", "tid|2|again");
			requireBadType(toIntList, server, silva, "tid|2|silva", "Life;Silva child");
			requireBadType(toIntList, server, ptHeader, "tid|2|pt");
			requireBadType(toIntList, server, scSilva, "tid|2|sc");
		}else{
			requireList(toIntList(toIntList, server, header, "tid|2|ncbi", "tid|2|again"), 2, 2);
			requireList(toIntList(toIntList, server, silva, "tid|2|silva", "Life;Silva child"), 2, 2);
			requireList(toIntList(toIntList, server, ptHeader, "tid|2|pt"), 2);
			requireList(toIntList(toIntList, server, scSilva, "tid|2|sc"), 2);
			requireMissing(toIntList, server, header, notFound, "tid|0|missing");
			requireMissing(toIntList, server, silva, notFound, "Unknown;Missing");
		}

		System.out.println("TAXSERVER_HEADER_ANCESTOR_TEST_PASS expect="
				+(baseline ? "baseline" : "fixed"));
	}

	/** Prevents construction of this standalone fixture. */
	private TaxServerHeaderAncestorTest(){}

	/** Builds the private TaxServer route-alias map without constructing a server.
	 * @return TaxServer.makeTypeMap() result
	 * @throws Exception If reflection cannot access or run makeTypeMap */
	@SuppressWarnings("unchecked")
	private static HashMap<String, Integer> typeMap() throws Exception{
		final Method method=TaxServer.class.getDeclaredMethod("makeTypeMap");
		method.setAccessible(true);
		return (HashMap<String, Integer>)method.invoke(null);
	}

	/** Allocates a TaxServer while bypassing its listener-starting constructor.
	 * @param tree Tiny in-memory tree
	 * @param notFound Counter for unresolved lookups
	 * @return TaxServer with only toIntList dependencies initialized
	 * @throws Exception If reflection cannot access Unsafe or private fields */
	private static TaxServer server(TaxTree tree, AtomicLong notFound) throws Exception{
		Class<?> unsafeClass=Class.forName("sun.misc.Unsafe");
		Field unsafeField=unsafeClass.getDeclaredField("theUnsafe");
		unsafeField.setAccessible(true);
		Object unsafe=unsafeField.get(null);
		Method allocateInstance=unsafeClass.getMethod("allocateInstance", Class.class);

		TaxServer server=(TaxServer)allocateInstance.invoke(unsafe, TaxServer.class);
		setField(server, "tree", tree);
		setField(server, "notFound", notFound);
		return server;
	}

	/** Sets one private TaxServer instance field.
	 * @param server Server instance to modify
	 * @param name Field name
	 * @param value Field value
	 * @throws Exception If the field cannot be reflected or set */
	private static void setField(TaxServer server, String name, Object value) throws Exception{
		final Field field=TaxServer.class.getDeclaredField(name);
		field.setAccessible(true);
		field.set(server, value);
	}

	/** Builds the smallest tree needed by old-style NCBI and SILVA headers.
	 * @return Tiny in-memory taxonomy tree with one resolvable child node */
	private static TaxTree tinyTree(){
		TaxNode[] nodes=new TaxNode[3];
		nodes[1]=new TaxNode(1, 1, TaxTree.LIFE, TaxTree.LIFE_E, "Life");
		nodes[2]=new TaxNode(2, 1, TaxTree.SPECIES, TaxTree.SPECIES_E, "Silva child");
		TaxTree tree=new TaxTree(nodes, null, 0, false, false, false, TaxTree.NO_RANK_E);
		tree.hashNames(false);
		return tree;
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
		return (IntList)method.invoke(server, type, (Object)names);
	}

	/** Requires an unsupported route to throw the exact baseline error.
	 * @param method Reflected toIntList method
	 * @param server Initialized server
	 * @param type Encoded route type
	 * @param names Query names
	 * @throws Exception If reflection fails */
	private static void requireBadType(Method method, TaxServer server, int type,
			String... names) throws Exception{
		try{
			toIntList(method, server, type, names);
		}catch(InvocationTargetException e){
			final Throwable cause=e.getCause();
			check(cause instanceof RuntimeException,
					"Expected RuntimeException; found "+cause);
			check(BAD_TYPE.equals(cause.getMessage()),
					"Expected "+BAD_TYPE+"; found "+cause.getMessage());
			return;
		}
		throw new AssertionError("Expected bad type for route "+type);
	}

	/** Requires a missing lookup to return no TaxIDs and increment notFound once.
	 * @param method Reflected toIntList method
	 * @param server Initialized server
	 * @param type Encoded route type
	 * @param notFound Missing-lookup counter
	 * @param name Missing query term
	 * @throws Exception If toIntList cannot be invoked or throws */
	private static void requireMissing(Method method, TaxServer server, int type,
			AtomicLong notFound, String name) throws Exception{
		final long before=notFound.get();
		requireList(toIntList(method, server, type, name));
		check(notFound.get()==before+1,
				"Expected one notFound increment; before="+before+", after="+notFound.get());
	}

	/** Requires missing TaxID lookup to be empty and reports the current counter delta.
	 * @param method Reflected toIntList method
	 * @param server Initialized server
	 * @param type Encoded TaxID route type
	 * @param notFound Missing-lookup counter
	 * @throws Exception If toIntList cannot be invoked or throws */
	private static void observeMissingTaxid(Method method, TaxServer server, int type,
			AtomicLong notFound) throws Exception{
		final long before=notFound.get();
		requireList(toIntList(method, server, type, "0"));
		System.out.println("TAXSERVER_MISSING_TAXID not_found_delta="
				+(notFound.get()-before));
	}

	/** Requires an alias to exist and encode to the expected route value.
	 * @param typeMap TaxServer alias map
	 * @param alias route alias to inspect
	 * @param expectedRaw expected encoded route value
	 * @return Actual route value */
	private static int requireAlias(final HashMap<String, Integer> typeMap,
			final String alias, final int expectedRaw){
		final Integer actual=typeMap.get(alias);
		check(actual!=null, "Missing route alias "+alias);
		check(actual.intValue()==expectedRaw,
				alias+" raw route should be "+expectedRaw+"; found "+actual);
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
	 * @param condition condition required for test success
	 * @param message diagnostic message for failure */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** Baseline unsupported-route exception message. */
	private static final String BAD_TYPE="{\"error\": \"Bad type\"}";
	/** TaxID route. */
	private static final int EXPECTED_TAXID=3;
	/** NCBI header route. */
	private static final int EXPECTED_HEADER=4;
	/** SILVA header route. */
	private static final int EXPECTED_SILVAHEADER=7;
	/** Plaintext NCBI header route. */
	private static final int EXPECTED_PT_HEADER=20;
	/** Semicolon SILVA header route. */
	private static final int EXPECTED_SC_SILVAHEADER=39;
}
