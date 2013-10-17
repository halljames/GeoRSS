using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;
using System.Xml.Serialization;


namespace GeoRSS
{
    [XmlRoot(ElementName = "rss")] // , Namespace = "http://www.w3.org/2005/Atom"
    public class Geo
    {
        [XmlAttribute("version")]
        public string version = "2.0";
        //public string title;
        //public string link;
        //public DateTime updated = DateTime.Now.ToUniversalTime();
        //public string id;

        [XmlElement("channel")]
        public List<Channel> ChannelList;

        public void AddChannel(Channel channel)
        {
            if (ChannelList == null)
            {
                ChannelList = new List<Channel>();
            }
            ChannelList.Add(channel);

        }
        public XmlDocument ToXML()
        {

            using (MemoryStream stream = new MemoryStream())
            {
                // xmlns:georss="http://www.georss.org/georss" 
                // xmlns:gml="http://www.opengis.net/gml"
                //xmlns:geo="http://www.w3.org/2003/01/geo/wgs84_pos#" 
                //xmlns:kml="http://www.opengis.net/kml/2.2" 
                //xmlns:dc="http://purl.org/dc/elements/1.1/"
                var ns = new XmlSerializerNamespaces();
                ns.Add("georss", "http://www.georss.org/georss");
                ns.Add("gml", "http://www.opengis.net/gml");
                ns.Add("geo", "http://www.w3.org/2003/01/geo/wgs84_pos#");
                ns.Add("kml", "http://www.opengis.net/kml/2.2");
                ns.Add("dc", "http://purl.org/dc/elements/1.1/");


                XmlSerializer s = new XmlSerializer(this.GetType());
                Console.WriteLine("Testing for type: {0}", this.GetType());
                s.Serialize(XmlWriter.Create(stream), this, ns);
                stream.Flush();
                stream.Seek(0, SeekOrigin.Begin);
                //object o = s.Deserialize(XmlReader.Create(stream));
                //Console.WriteLine("  Deserialized type: {0}", o.GetType());
                XmlDocument xml = new XmlDocument();
                xml.Load(stream);
                Console.Write(xml.InnerXml);
                return xml;
            }

            //var serializer = new XmlSerializer(this.GetType());
            //serializer.Serialize(new StreamWriter("test.xml"), this);



        }
    }
    public class Channel
    {

        public string title;
        public string description;
        public string link;


        [XmlElement("item")]
        public List<Item> ItemList;


        public void AddItem(Item item)
        {
            if (ItemList == null)
            {
                ItemList = new List<Item>();
            }
            ItemList.Add(item);

        }
    }

    public class Item
    {
        public string title;
        public string description;

        //public string link;
        //public string id;


        public string pubDate;

        [XmlIgnore]
        public DateTime Date
        {
            get { return Date; }
            set { pubDate = value.ToString("r"); }
        }


        [XmlElement(ElementName = "point", Namespace = "http://www.georss.org/georss")]
        public string point;


        //[XmlIgnore]
        //public Double lat;
        //[XmlIgnore]
        //public Double lng;

        public void AddPoint(Double lat, Double lng)
        {
            this.point = Math.Round(lat, 3).ToString() + " " + Math.Round(lng, 3).ToString();
        }
        public void AddPointEastingNorthing(string easting, string northing)
        {
            GeoUtils util = new GeoUtils();
            util.ne2ll(easting, northing,1,6);
            AddPoint(Math.Round(Convert.ToDouble(util.LatitudeValue), 4), Math.Round(Convert.ToDouble(util.LongtitudeValue), 4));

        }

    }
}
