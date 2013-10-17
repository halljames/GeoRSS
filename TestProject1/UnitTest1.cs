using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using GeoRSS;
using System.Xml;
namespace TestProject1
{
    [TestClass]
    public class UnitTest1
    {
        [TestMethod]
        public void TestMethod1()
        {
            Geo test = new Geo();


            Channel channel = new Channel();
            channel.title = "Test RSS";
            channel.description = "Test Geo RSS";
            channel.link = "http://www.google.co.uk";

            Item item = new Item();
            item.title = "Test Entry 1";
            item.description = "Bla Bla Bla";
            item.Date = DateTime.Now;
            // 54.898742,-1.398246
            item.AddPoint(54.898742, -1.398246);
            //item.AddPoint(434171, 550112);
            channel.AddItem(item);
            item = new Item();
            item.title = "Test Entry 2";
            item.description = "Lla Lla Lla";
            item.Date = DateTime.Now;
            item.AddPoint(53.456,-1.345);
            //item.AddPoint(434171, 530112);
            channel.AddItem(item);



            test.AddChannel(channel);



            XmlDocument xml = test.ToXML();
            xml.Save("C:\\Temp\\rss.xml");
        }
    }
}
