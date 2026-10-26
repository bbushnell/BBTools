package test.map;

import java.lang.reflect.Method;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import map.LongArrayListHashMap;
import map.LongLongListHashMap;
import structures.LongList;

/** Ownership and content checks for copying insertion and no-copy list-map growth.
 * @author Shinobu
 */
public final class ListMapResizeTest{

	/** Captures old alias detachment or verifies repaired ownership and adoption.
	 * @param args One expect=baseline|fixed argument
	 * @throws Exception If reflective access to the new public API fails */
	public static void main(String[] args) throws Exception{
		check(args.length==1, "Expected expect=baseline|fixed");
		final boolean fixed=args[0].equals("expect=fixed");
		check(fixed || args[0].equals("expect=baseline"), "Invalid expectation");
		array(fixed);
		longs(fixed);
		System.out.println("PASS list-map ownership fixed="+fixed);
	}

	/** Checks object-list copying, growth aliases, replacement and recycled collision slots. */
	private static void array(final boolean fixed) throws Exception{
		final LongArrayListHashMap<Long> map=new LongArrayListHashMap<Long>(4, 0.7f);
		final ArrayList<Long> input=arrayList(7, 8);
		map.put(1, input);
		final ArrayList<Long> owned=map.get(1);
		final ArrayList<Long>[] oldArray=map.values();
		check(owned!=input, "Ordinary put must copy caller-owned object list");
		input.clear();
		final ArrayList<Long> extra=arrayList(3, 4);
		map.put(1, extra);
		extra.clear();
		check(owned.size()==4 && owned.get(2)==3L && owned.get(3)==4L, "Existing-key bulk put copies and appends");
		for(int i=2; i<128; i++){map.put(i, (long)i);}
		check(map.values()!=oldArray, "Backing value array must grow");
		check((map.get(1)==owned)==fixed, "Object-list identity across resize");
		map.put(1, 9L);
		check(map.get(1).size()==5 && map.get(1).get(0)==7L && map.get(1).get(4)==9L, "Object-list contents after resize");
		check(owned.size()==(fixed ? 5 : 4), "Earlier alias observes append only when retained");
		check(map.verify(), "Healthy object-list map");
		//Brian retained empty-list rejection in verify(), not insertion. Growth must still preserve an empty input.
		map.put(500, new ArrayList<Long>());
		for(int i=128; i<1024; i++){if(i!=500){map.put(i, (long)i);}}
		check(map.get(500)!=null && map.get(500).isEmpty(), "Empty object-list mapping survives growth");
		Method adopt=null;
		try{adopt=LongArrayListHashMap.class.getMethod("putCopyFree", long.class, ArrayList.class);}
		catch(NoSuchMethodException e){check(!fixed, "Missing object-list adoption method");}
		if(!fixed){check(adopt==null, "Baseline must lack adoption API"); return;}
		final LongArrayListHashMap<Long> collision=new LongArrayListHashMap<Long>(4, 0.7f);
		final long step=collision.keys().length-10;
		final ArrayList<Long> first=arrayList(10, 11), second=arrayList(20, 21), third=arrayList(30, 31);
		check((Boolean)adopt.invoke(collision, 0L, first), "New object key");
		adopt.invoke(collision, step, second);
		adopt.invoke(collision, 2*step, third);
		check(collision.get(0)==first && collision.get(step)==second, "Exact object lists adopted");
		final ArrayList<Long> replacement=arrayList(40, 41);
		check(!(Boolean)adopt.invoke(collision, step, replacement) && collision.size()==3, "Replacement is not append or insertion");
		check(collision.get(step)==replacement && second.size()==2, "Detached old object list remains unchanged");
		check(!(Boolean)adopt.invoke(collision, step, replacement), "Same-reference replacement");
		check(collision.remove(0) && first.isEmpty(), "Removal clears adopted object list");
		check(collision.get(step)==replacement && collision.get(2*step)==third, "Collision rehash retains surviving aliases");
		collision.put(3*step, 50L);
		check(collision.get(3*step)==first && first.get(0)==50L, "Removed small list may be recycled for another key");
		check(collision.verify(), "Object-list collision invariants");
		collision.clear();
		check(collision.isEmpty() && first.isEmpty() && replacement.isEmpty() && third.isEmpty(), "Clear owns adopted lists");
		final ArrayList<Long> empty=new ArrayList<Long>();
		check((Boolean)adopt.invoke(collision, 1000L, empty), "Direct empty object-list adoption");
		checkNullRejected(adopt, collision, 1001L);
		check(!collision.contains(1001), "Null rejection must not insert an object-list key");
		for(long i=2000; i<2032; i++){collision.put(i, Long.valueOf(i));}
		check(collision.get(1000)==empty && empty.isEmpty(), "Direct empty object-list adoption survives growth");
	}

	/** Checks the same contracts for primitive long lists. */
	private static void longs(final boolean fixed) throws Exception{
		final LongLongListHashMap map=new LongLongListHashMap(4, 0.7f);
		final LongList input=longList(7, 8);
		map.put(1, input);
		final LongList owned=map.get(1);
		final LongList[] oldArray=map.values();
		check(owned!=input, "Ordinary put must copy caller-owned primitive list");
		input.clear();
		final LongList extra=longList(3, 4);
		map.put(1, extra);
		extra.clear();
		check(owned.size==4 && owned.get(2)==3 && owned.get(3)==4, "Existing-key primitive bulk put copies and appends");
		for(int i=2; i<128; i++){map.put(i, (long)i);}
		check(map.values()!=oldArray, "Primitive backing array must grow");
		check((map.get(1)==owned)==fixed, "Primitive-list identity across resize");
		map.put(1, 9L);
		check(map.get(1).size==5 && map.get(1).get(0)==7 && map.get(1).get(4)==9, "Primitive-list contents after resize");
		check(owned.size==(fixed ? 5 : 4), "Earlier primitive alias observes append only when retained");
		check(map.verify(), "Healthy primitive-list map");
		map.put(500, new LongList(1));
		for(int i=128; i<1024; i++){if(i!=500){map.put(i, (long)i);}}
		check(map.get(500)!=null && map.get(500).isEmpty(), "Empty primitive-list mapping survives growth");
		Method adopt=null;
		try{adopt=LongLongListHashMap.class.getMethod("putCopyFree", long.class, LongList.class);}
		catch(NoSuchMethodException e){check(!fixed, "Missing primitive-list adoption method");}
		if(!fixed){check(adopt==null, "Baseline must lack adoption API"); return;}
		final LongLongListHashMap collision=new LongLongListHashMap(4, 0.7f);
		final long step=collision.keys().length-10;
		final LongList first=longList(10, 11), second=longList(20, 21), third=longList(30, 31);
		check((Boolean)adopt.invoke(collision, 0L, first), "New primitive key");
		adopt.invoke(collision, step, second);
		adopt.invoke(collision, 2*step, third);
		check(collision.get(0)==first && collision.get(step)==second, "Exact primitive lists adopted");
		final LongList replacement=longList(40, 41);
		check(!(Boolean)adopt.invoke(collision, step, replacement) && collision.size()==3, "Primitive replacement does not append");
		check(collision.get(step)==replacement && second.size==2, "Detached old primitive list remains unchanged");
		check(!(Boolean)adopt.invoke(collision, step, replacement), "Same-reference primitive replacement");
		check(collision.remove(0) && first.isEmpty(), "Removal clears adopted primitive list");
		check(collision.get(step)==replacement && collision.get(2*step)==third, "Primitive rehash retains aliases");
		collision.put(3*step, 50L);
		check(collision.get(3*step)==first && first.get(0)==50, "Primitive recycled slot");
		check(collision.verify(), "Primitive collision invariants");
		collision.clear();
		check(collision.isEmpty() && first.isEmpty() && replacement.isEmpty() && third.isEmpty(), "Clear owns primitive lists");
		final LongList empty=new LongList(1);
		check((Boolean)adopt.invoke(collision, 1000L, empty), "Direct empty primitive-list adoption");
		checkNullRejected(adopt, collision, 1001L);
		check(!collision.contains(1001), "Null rejection must not insert a primitive-list key");
		for(long i=2000; i<2032; i++){collision.put(i, i);}
		check(collision.get(1000)==empty && empty.isEmpty(), "Direct empty primitive-list adoption survives growth");
	}

	/** Invalid null input must fail before creating an occupied mapping, even under -da. */
	private static void checkNullRejected(final Method adopt, final Object map, final long key) throws Exception{
		boolean failed=false;
		try{adopt.invoke(map, key, null);}catch(InvocationTargetException e){failed=e.getCause() instanceof NullPointerException;}
		check(failed, "Adoption must reject null before mutation");
	}

	/** Creates an independent two-value object list. */
	private static ArrayList<Long> arrayList(final long a, final long b){
		final ArrayList<Long> list=new ArrayList<Long>(2);
		list.add(a); list.add(b);
		return list;
	}

	/** Creates an independent two-value primitive list. */
	private static LongList longList(final long a, final long b){
		final LongList list=new LongList(2);
		list.add(a); list.add(b);
		return list;
	}

	/** Keeps checks active with assertions disabled. */
	private static void check(final boolean ok, final String message){
		if(!ok){throw new AssertionError(message);}
	}

	/** Static fixture only. */
	private ListMapResizeTest(){}
}
