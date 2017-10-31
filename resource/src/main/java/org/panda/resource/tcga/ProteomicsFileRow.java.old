package org.panda.resource.tcga;

import org.panda.resource.PhosphoSitePlus;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * This class represents a single row in a proteimics data file. This also includes TCGA RPPA files.
 *
 * @author Ozgun Babur
 */
public class ProteomicsFileRow implements Cloneable
{
	public String id;
	public double[] vals;
	public List<String> genes;
	public Map<String, List<String>> sites;
	public SiteEffect effect;
	public String[] header;
	Type type;

	public ProteomicsFileRow(String id, double[] vals, List<String> genes, Map<String, List<String>> sites)
	{
		this.id = id;
		this.vals = vals;
		this.genes = genes;
		this.sites = sites;

		if (sites != null && !sites.isEmpty()) type = Type.SITE_SPECIFIC;
		else type = Type.TOTAL_PROTEIN;

		if (sites == null) return;

		for (String gene : sites.keySet())
		{
			for (String site : sites.get(gene))
			{
				Integer eff = PhosphoSitePlus.get().getEffect(gene, site);
				if (eff != null)
				{
					effect = SiteEffect.getValue(eff);
					break;
				}
			}
			if (effect != null) break;
		}

//		if (effect == null && !sites.isEmpty())
//		{
//			for (String site : sites)
//			{
//				for (String gene : genes)
//				{
//					Integer eff = PhosphoSitePlus.getClosestEffect(gene, site);
//					if (eff != null)
//					{
//						effect = SiteEffect.getValue(eff);
//						break;
//					}
//				}
//			}
//		}
	}

	public boolean isPhospho()
	{
		return type == Type.SITE_SPECIFIC;
	}

	public boolean isTotalProt()
	{
		return type == Type.TOTAL_PROTEIN;
	}

	public boolean isActivity()
	{
		return type == Type.ACTIVITY;
	}

	public int getSelfEffect()
	{
		if (!isPhospho()) return 1;
		else if (effect == null) return 0;
		else return effect.getVal();
	}

	@Override
	public boolean equals(Object obj)
	{
		return obj instanceof ProteomicsFileRow && id.equals(((ProteomicsFileRow) obj).id);
	}

	@Override
	public int hashCode()
	{
		return id.hashCode();
	}

	public enum SiteEffect
	{
		ACTIVATING(1),
		INHIBITING(-1),
		COMPLEX(0);

		int val;

		private SiteEffect(int val)
		{
			this.val = val;
		}

		public int getVal()
		{
			return val;
		}

		public static SiteEffect getValue(int x)
		{
			for (SiteEffect effect : values())
			{
				if (effect.val == x) return effect;
			}
			return null;
		}
	}

	public void makeActivityNode(boolean isActivated)
	{
		type = Type.ACTIVITY;
		vals = new double[1];
		vals[0] = isActivated ? 1 : -1;
		effect = null;
		sites = null;
	}

	@Override
	public Object clone()
	{
		try
		{
			ProteomicsFileRow clone = (ProteomicsFileRow) super.clone();

//			clone.vals = vals.clone();

			return clone;
		}
		catch (CloneNotSupportedException e)
		{
			throw new RuntimeException(e);
		}
	}

	@Override
	public String toString()
	{
		return id;
	}


	public static void write(Collection<ProteomicsFileRow> datas, String filename)
	{
		try
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

			writer.write("ID\tSymbol\tSite\tEffect");

			Iterator<ProteomicsFileRow> iter = datas.iterator();
			ProteomicsFileRow sample = iter.next();
			while (sample.isActivity()) sample = iter.next();

			for (int i = 0; i < sample.vals.length; i++)
			{
				writer.write(sample.header == null ? "\tv" + i : "\t" + sample.header[i]);
			}

			for (ProteomicsFileRow data : datas)
			{
				if (data.isActivity()) continue;

				writer.write("\n" + data.id);

				String items = new ArrayList<>(data.genes).toString().
					replace(",", "").replace("[", "").replace("]", "");
				writer.write("\t" + items);

				if (data.sites != null)
				{
					List<String> siteList = new ArrayList<String>();
					for (String gene : data.sites.keySet())
					{
						List<String> ss = data.sites.get(gene);
						String s = ss.toString().replace("[", "").replace("]", "").replace(", ", "|");
						siteList.add(s);
					}

					items = siteList.toString().replace(",", "").replace("[", "").replace("]", "");
					writer.write("\t" + items);

					writer.write("\t" + (data.effect == null ? "" : data.effect == SiteEffect.COMPLEX ?
						"c" : data.effect == SiteEffect.ACTIVATING ? "a" : "i"));
				}
				else writer.write("\t\t");

				for (double v : data.vals)
				{
					writer.write("\t" + v);
				}
			}

			writer.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public enum Type
	{
		TOTAL_PROTEIN,
		SITE_SPECIFIC,
		ACTIVITY,
		EXPRESSION
	}

	public static List<ProteomicsFileRow> copy(List<ProteomicsFileRow> orig)
	{
		List<ProteomicsFileRow> list = new ArrayList<ProteomicsFileRow>();
		for (ProteomicsFileRow data : orig)
		{
			list.add((ProteomicsFileRow) data.clone());
		}
		return list;
	}
}
