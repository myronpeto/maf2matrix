package org.panda.resource;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Provides genes in OncoKB.
 *
 * How to update this resource:
 * Go to oncokb.org and click on the "genes" link. All genes should come in one page. Copy everything and paste into
 * a text-only editor. Delete unnecessary portion.
 *
 * @author Ozgun Babur
 */
public class OncoKB extends FileServer
{
	private static OncoKB instance;

	private Map<String, String> sym2level;

	public static synchronized OncoKB get()
	{
		if (instance == null) instance = new OncoKB();
		return instance;
	}

	public Set<String> getAllSymbols()
	{
		return sym2level.keySet();
	}

	public boolean isCancerGene(String sym)
	{
		return sym2level.containsKey(sym);
	}

	public String getLevel(String gene)
	{
		return sym2level.get(gene);
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"oncokb.txt"};
	}

	@Override
	public boolean load() throws IOException
	{
		sym2level = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0]).skip(1).map(l -> l.split("\t")).forEach(t ->
				sym2level.put(t[0], t[1]));

		return true;
	}

	public static void main(String[] args)
	{
		System.out.println(get().isCancerGene("AR"));
	}
}
